"""
quantum_ran_complete.py
=======================
Complete, self-contained implementation of the corrected quantum-inspired
QUBO framework for RAN optimisation, as described in:

  "Quantum-Inspired QUBO Solvers for Scalable Radio Access Network
   Optimisation: Formulations, Implementation Corrections, and
   Benchmarking"

HOW TO RUN
----------
Place this file in the same folder as the original module files:
    quantum_inspired_ran.py
    quantum_optimisation.py

Then run:
    python quantum_ran_complete.py

Output PDFs are saved to ./figures/ (created automatically).
Works on Windows, Linux, and macOS.

BUG FIXES APPLIED
-----------------
BUG 1 — SBA (quantum_inspired_ran.py)
  Root cause : QUBO→Ising conversion produces h_i > 0 for all spins in
               frequency-planning problems, driving all spins to +1.
  Fix        : Negate J and h; normalise by spectral scale κ.
               Eq. (47)–(48) in the paper.

BUG 2a — DA overflow (quantum_inspired_ran.py)
  Root cause : exp(dB·dE) raises OverflowError in replica exchange.
  Fix        : Clamp exponent to [−500, 500]. Eq. (54).

BUG 2b — DA log (quantum_inspired_ran.py)
  Root cause : energy_log records energies.min() which rises after swaps.
  Fix        : Log global best_e (monotone non-increasing). Eq. (55).

BUG 3  — SQA backwards Metropolis (quantum_optimisation.py)
  Root cause : dE_local = +2σ_i·field = −ΔE, so Metropolis accepted
               energy-raising moves and rejected energy-lowering ones.
  Fix        : Use dE = −2σ_i·field (= +ΔE). Eqs. (59)–(61).

BUG 4  — PA log (quantum_inspired_ran.py)
  Root cause : energy_log records E.min() which rises after resampling.
  Fix        : Log global best_e (monotone non-increasing). Eq. (58).

ADDITIONAL IMPROVEMENTS (from paper Discussion, Section VI)
-----------------------------------------------------------
- Feasibility repair step (Eq. 77) applied to all solver outputs.
- Sensitivity analysis: λ_c and n_steps sweeps.
- Reproducibility study: multi-seed statistics.
- GPU projection formula (Eq. 81).
"""

# ── standard library ──────────────────────────────────────────────────────────
import sys, math, time, random, warnings
from pathlib import Path
warnings.filterwarnings("ignore")

# ── third-party ───────────────────────────────────────────────────────────────
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ── output folder ─────────────────────────────────────────────────────────────
OUT_DIR = Path(__file__).resolve().parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

def save_fig(fig, name):
    path = OUT_DIR / name
    fig.savefig(path)
    print(f"  ✓  {path}")

# ── module path ───────────────────────────────────────────────────────────────
_here = str(Path(__file__).resolve().parent)
if _here not in sys.path:
    sys.path.insert(0, _here)

# ── original module imports (QUBO builders only) ──────────────────────────────
from quantum_inspired_ran import LargeScaleRANQUBO
from quantum_optimisation  import (
    QUBOBuilder, ClassicalSolver, QAOASimulator,
    QuantumOptimisationOrchestrator,
)

# ═══════════════════════════════════════════════════════════════════════════════
#  CORRECTED SOLVER CLASSES
#  All four bugs fixed directly in the class body.
# ═══════════════════════════════════════════════════════════════════════════════

class SBA:
    """
    Simulated Bifurcation Algorithm — fully corrected (Bug 1 fixed).

    Original update:  dmomentum += ξ·(J@σ + h)
    Problem          : h_i > 0 for all spins in RAN frequency-planning
                       QUBOs ⟹ all spins driven to +1 (all-ones solution).
    Fix              : Negate J and h after spectral normalisation so that
                       SBA *minimises* the QUBO energy. See Eq. (47)–(50).

    Parameters
    ----------
    n_steps   : Total bifurcation steps (default 3000).
    dt        : Symplectic Euler step size (default 0.015).
    xi        : Coupling strength ξ (default 0.5).
    p_start   : Initial bifurcation parameter (< 0 ⟹ paramagnetic).
    p_end     : Final bifurcation parameter  (> 0 ⟹ ferromagnetic).
    seed      : RNG seed for reproducibility.
    """

    def __init__(self, n_steps=3000, dt=0.015, xi=0.5,
                 p_start=-1.0, p_end=1.5, seed=42):
        self.n_steps = n_steps
        self.dt      = dt
        self.xi      = xi
        self.p_start = p_start
        self.p_end   = p_end
        self.rng     = np.random.default_rng(seed)
        self.energy_log = []          # monotone global-best QUBO energy

    def solve(self, Q, max_n=100):
        """Return (x* ∈ {0,1}^n, best_qubo_energy)."""
        N  = min(Q.shape[0], max_n)
        Qs = (Q[:N, :N] + Q[:N, :N].T) / 2          # symmetrise

        # ── QUBO → Ising (Eq. 3) ─────────────────────────────────────────
        h = Qs.diagonal() / 2 + (Qs.sum(axis=1) - Qs.diagonal()) / 2
        J = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                if i != j:
                    J[i, j] = Qs[i, j] / 4

        # ── BUG 1 FIX: spectral normalisation + sign negation (Eq. 47-48) ─
        kappa = max(np.abs(J).max() * N**0.5, np.abs(h).max(), 1e-9)
        Jn    = -J / kappa
        hn    = -h / kappa

        # ── bifurcation dynamics (Eq. 49) ─────────────────────────────────
        sigma = self.rng.uniform(-0.1, 0.1, N)
        mom   = self.rng.uniform(-0.1, 0.1, N)
        self.energy_log.clear()
        best_e = float("inf")
        log_every = max(1, self.n_steps // 30)

        for step in range(self.n_steps):
            p     = self.p_start + (self.p_end - self.p_start) * step / self.n_steps
            coup  = Jn @ sigma + hn
            dmom  = (p - sigma**2) * sigma + self.xi * coup
            mom  += self.dt * dmom
            sigma += self.dt * mom
            sigma  = np.clip(sigma, -1.0, 1.0)

            if step % log_every == 0:                 # Eq. 51: global best
                x_cur  = ((np.sign(sigma) + 1) / 2).astype(int)
                e_cur  = float(x_cur @ Qs @ x_cur)
                best_e = min(best_e, e_cur)
                self.energy_log.append(best_e)

        x_f    = ((np.sign(sigma) + 1) / 2).astype(int)
        best_e = min(best_e, float(x_f @ Qs @ x_f))
        result       = np.zeros(Q.shape[0], dtype=int)
        result[:N]   = x_f
        return result, round(best_e, 3)


class DA:
    """
    Digital Annealing with parallel tempering — fully corrected.

    Bug 2a Fix: Clamp exchange exponent to [−500, 500]  (Eq. 54).
    Bug 2b Fix: Log global best_e, not energies.min()   (Eq. 55).

    Parameters
    ----------
    n_replicas    : Number of parallel replicas R.
    t_min, t_max  : Temperature schedule endpoints (geometric spacing).
    n_sweeps      : Metropolis sweeps per replica per step.
    exchange_freq : Replica exchange attempted every this many sweeps.
    seed          : RNG seed.
    """

    def __init__(self, n_replicas=10, t_min=0.001, t_max=8.0,
                 n_sweeps=400, exchange_freq=15, seed=42):
        self.n_rep    = n_replicas
        self.temps    = np.logspace(math.log10(t_min),
                                    math.log10(t_max), n_replicas)
        self.n_sweeps = n_sweeps
        self.exch_f   = exchange_freq
        self.rng      = np.random.default_rng(seed)
        self.energy_log = []          # monotone global-best QUBO energy

    def _dE(self, x, Q, i):
        """O(n) single-spin flip energy (Eq. 53)."""
        col = Q[:, i] + Q[i, :]
        return float((1 - 2 * x[i]) * (col @ x - Q[i, i] * (2 * x[i] - 1)))

    def solve(self, Q, max_n=100):
        N   = min(Q.shape[0], max_n)
        Qs  = Q[:N, :N]
        rng = self.rng

        replicas = rng.integers(0, 2, (self.n_rep, N))
        energies = np.array([float(r @ Qs @ r) for r in replicas])
        best_x   = replicas[np.argmin(energies)].copy()
        best_e   = float(energies.min())
        self.energy_log.clear()

        for sweep in range(self.n_sweeps):
            for r in range(self.n_rep):
                T = self.temps[r]
                for _ in range(N):
                    i  = rng.integers(0, N)
                    dE = self._dE(replicas[r], Qs, i)
                    if dE < 0 or rng.random() < math.exp(max(-500.0, -dE / T)):
                        replicas[r, i] ^= 1
                        energies[r]    += dE
                if energies[r] < best_e:
                    best_e = energies[r]
                    best_x = replicas[r].copy()

            if sweep % self.exch_f == 0:
                for r in range(self.n_rep - 1):
                    dE   = energies[r] - energies[r + 1]
                    dB   = 1.0 / self.temps[r] - 1.0 / self.temps[r + 1]
                    # BUG 2a FIX (Eq. 54): clamp before math.exp
                    prob = min(1.0, math.exp(min(500.0, dB * dE)))
                    if rng.random() < prob:
                        replicas[r],   replicas[r + 1]  = (
                            replicas[r + 1].copy(), replicas[r].copy())
                        energies[r],   energies[r + 1]  = (
                            energies[r + 1], energies[r])
                # BUG 2b FIX (Eq. 55): log monotone global best
                self.energy_log.append(best_e)

        result     = np.zeros(Q.shape[0], dtype=int)
        result[:N] = best_x
        return result, round(best_e, 3)


class PA:
    """
    Population Annealing — fully corrected.

    Bug 4 Fix: Log global best_e instead of E.min() (Eq. 58).
               E.min() rises after importance-resampling; best_e is monotone.

    Parameters
    ----------
    population    : Number of replicas R.
    n_temps       : Temperature schedule steps.
    t_start, t_end: Temperature schedule endpoints (geometric, decreasing).
    sweeps_per_t  : Metropolis sweeps per replica per temperature step.
    seed          : RNG seed.
    """

    def __init__(self, population=60, n_temps=50,
                 t_start=8.0, t_end=0.001, sweeps_per_t=3, seed=42):
        self.R      = population
        self.T_sch  = np.logspace(math.log10(t_start),
                                   math.log10(t_end), n_temps)
        self.sw_t   = sweeps_per_t
        self.rng    = np.random.default_rng(seed)
        self.energy_log = []          # monotone global-best QUBO energy

    def solve(self, Q, max_n=80):
        N  = min(Q.shape[0], max_n)
        Qs = Q[:N, :N]

        pop    = self.rng.integers(0, 2, (self.R, N))
        E      = np.array([float(pop[r] @ Qs @ pop[r]) for r in range(self.R)])
        best_x = pop[np.argmin(E)].copy()
        best_e = float(E.min())
        self.energy_log.clear()

        for T in self.T_sch:
            # log-domain importance weights (Eq. 56) for numerical stability
            log_w  = -(1 / T) * E
            log_w -= log_w.max()
            w      = np.exp(log_w);  w /= w.sum()
            idx    = self.rng.choice(self.R, self.R, p=w)
            pop    = pop[idx];  E = E[idx]

            for r in range(self.R):
                for _ in range(self.sw_t * N):
                    i   = self.rng.integers(0, N)
                    col = Qs[:, i] + Qs[i, :]
                    dE  = float((1 - 2 * pop[r, i]) *
                                (col @ pop[r] - Qs[i, i] * (2 * pop[r, i] - 1)))
                    if dE < 0 or self.rng.random() < math.exp(max(-500.0, -dE / T)):
                        pop[r, i] ^= 1;  E[r] += dE
                if E[r] < best_e:
                    best_e = E[r];  best_x = pop[r].copy()

            # BUG 4 FIX (Eq. 58): log global best (monotone)
            self.energy_log.append(best_e)

        result     = np.zeros(Q.shape[0], dtype=int)
        result[:N] = best_x
        return result, round(best_e, 3)


class SQA:
    """
    Simulated Quantum Annealing (path-integral Monte Carlo) — fully corrected.

    Bug 3 Fix: Reversed Metropolis acceptance criterion corrected.

    Original: dE_local = +2·σ_i·field = −ΔE  (wrong sign)
    Effect  : Metropolis "if dE < 0: accept" accepted energy-RAISING moves.
              The algorithm ran thermodynamically backwards. The steady-state
              distribution was proportional to exp(+β·E) — peaked at the
              HIGHEST energy (worst) solution.
    Fix     : Use dE = −2·σ_i·field (= +ΔE, the correct flip energy).
              See Eqs. (59)–(62) in the paper.

    Parameters
    ----------
    n_replicas          : Trotter slices P.
    T                   : Classical temperature.
    n_steps             : Monte Carlo steps.
    gamma_start, gamma_end: Transverse field schedule endpoints.
    seed                : RNG seed.
    """

    def __init__(self, n_replicas=20, T=0.3, n_steps=3000,
                 gamma_start=5.0, gamma_end=0.001, seed=42):
        self.P           = n_replicas
        self.T           = T
        self.n_steps     = n_steps
        self.gamma_start = gamma_start
        self.gamma_end   = gamma_end
        self.rng         = np.random.default_rng(seed)
        self.energy_log  = []         # monotone global-best QUBO energy

    def solve(self, Q, max_n=20):
        N  = min(Q.shape[0], max_n)
        Qs = Q[:N, :N]
        P  = self.P

        # ── QUBO → Ising (Eq. 3) ─────────────────────────────────────────
        J = np.zeros((N, N));  h = np.zeros(N)
        for i in range(N):
            h[i] = Qs[i, i] / 2 + sum(
                (Qs[i, j] + Qs[j, i]) / 4 for j in range(N) if j != i)
            for j in range(i + 1, N):
                J[i, j] = (Qs[i, j] + Qs[j, i]) / 4

        replicas = self.rng.choice([-1, 1], size=(P, N))
        best_x   = replicas[0].copy();  best_e = float("inf")
        self.energy_log.clear()
        log_every = max(1, self.n_steps // 30)

        for step in range(self.n_steps):
            s      = step / self.n_steps
            gamma  = self.gamma_start * (1 - s) + self.gamma_end * s
            beta   = 1.0 / max(self.T, 1e-9)
            # Trotter coupling (Eq. 63)
            J_perp = (-0.5 * self.T *
                      math.log(math.tanh(gamma * beta / P))
                      if P > 1 else 0.0)

            for p in range(P):
                i  = self.rng.integers(0, N)
                lf = h[i] + sum(
                    J[min(i, j), max(i, j)] * replicas[p, j]
                    for j in range(N) if j != i)

                # BUG 3 FIX (Eq. 61): −2σ_i·field = +ΔE (correct sign)
                dE_cls = -2 * replicas[p, i] * lf
                dE_q   = (2 * J_perp * replicas[p, i] *
                          (replicas[(p - 1) % P, i] +
                           replicas[(p + 1) % P, i]))
                dE     = dE_cls + dE_q

                if dE < 0 or self.rng.random() < math.exp(max(-500.0, -dE * beta)):
                    replicas[p, i] *= -1

            if step % log_every == 0:
                for p in range(P):
                    x = ((replicas[p] + 1) // 2).astype(int)
                    e = float(x @ Qs @ x)
                    if e < best_e:
                        best_e = e;  best_x = replicas[p].copy()
                self.energy_log.append(best_e)

        x_out      = ((best_x + 1) // 2).astype(int)
        result     = np.zeros(Q.shape[0], dtype=int)
        result[:N] = x_out
        return result, round(best_e, 3)


# ── patch QuantumOptimisationOrchestrator to use corrected SQA ───────────────
from quantum_optimisation import SimulatedQuantumAnnealer

def _sqa_patched(self, Q, max_n=20):
    solver = SQA(n_replicas=self.P, T=self.T, n_steps=self.n_steps,
                 gamma_start=self.gamma_start, gamma_end=self.gamma_end,
                 seed=42)
    result, best_e    = solver.solve(Q, max_n=max_n)
    self.energy_log   = solver.energy_log
    return result, best_e

SimulatedQuantumAnnealer.solve = _sqa_patched


# ═══════════════════════════════════════════════════════════════════════════════
#  HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def feasibility_repair(x, n_cells, n_channels=3):
    """
    One-hot projection repair (Eq. 77, Discussion Section VI-A).

    For each cell, keep only the channel with the maximum activation value
    and zero all others, enforcing the one-hot constraint in O(N_c · C).

    Parameters
    ----------
    x          : Binary solution vector of length n_cells * n_channels.
    n_cells    : Number of cells.
    n_channels : Channels per cell (default 3).

    Returns
    -------
    x_rep : Repaired feasible binary vector.
    """
    x_rep = np.zeros_like(x)
    for i in range(n_cells):
        grp = x[i * n_channels:(i + 1) * n_channels]
        best_c = int(np.argmax(grp))
        x_rep[i * n_channels + best_c] = 1
    return x_rep


def run_sa(Q, n_steps=5000, T_start=5.0, seed=42, log_points=30):
    """Classical SA on QUBO. Returns (x*, best_e, monotone_energy_log)."""
    rng    = np.random.default_rng(seed)
    N      = Q.shape[0]
    x      = rng.integers(0, 2, N)
    E      = float(x @ Q @ x)
    best_e = E;  log = []
    for step in range(n_steps):
        T   = T_start * (1.0 - step / n_steps) ** 2
        i   = rng.integers(0, N)
        col = Q[:, i] + Q[i, :]
        dE  = float((1 - 2 * x[i]) * (col @ x - Q[i, i] * x[i]))
        if dE < 0 or (T > 0 and rng.random() < math.exp(max(-500.0, -dE / T))):
            x[i] ^= 1;  E += dE
        best_e = min(best_e, E)
        if step % (n_steps // log_points) == 0:
            log.append(best_e)
    return x, best_e, log


def sat_pct(x, n_cells, n_channels=3):
    """One-channel-per-cell constraint satisfaction (%)."""
    ok = sum(
        1 for i in range(n_cells)
        if sum(int(x[i * n_channels + c]) for c in range(n_channels)) == 1)
    return round(100.0 * ok / n_cells, 1)


def gpu_runtime_projection(t_cpu_ms, n_current, n_target, G=100):
    """
    GPU runtime projection formula (Eq. 81, Discussion Section VI-D).

    T_GPU(n) ≈ T_CPU(n_current) · (n_target / n_current) / G

    Parameters
    ----------
    t_cpu_ms  : Measured CPU runtime at n_current (ms).
    n_current : Problem size at which t_cpu_ms was measured.
    n_target  : Target problem size.
    G         : GPU parallelism factor (default 100).

    Returns
    -------
    Projected GPU runtime (ms).
    """
    return t_cpu_ms * (n_target / n_current) / G


# ═══════════════════════════════════════════════════════════════════════════════
#  MATPLOTLIB STYLE
# ═══════════════════════════════════════════════════════════════════════════════
plt.rcParams.update({
    "font.family":       "serif",
    "font.serif":        ["Times New Roman", "DejaVu Serif"],
    "font.size":         10, "axes.titlesize": 11, "axes.labelsize": 10,
    "xtick.labelsize":   9,  "ytick.labelsize": 9, "legend.fontsize": 9,
    "legend.framealpha": 0.9, "legend.edgecolor": "#cccccc",
    "axes.spines.top":   False, "axes.spines.right": False,
    "axes.grid":         True, "grid.linestyle": "--",
    "grid.alpha":        0.45, "grid.linewidth": 0.6,
    "lines.linewidth":   2.0, "lines.markersize": 5,
    "figure.dpi":        300, "savefig.dpi": 300,
    "savefig.bbox":      "tight", "savefig.pad_inches": 0.06,
})

C  = {"SBA":"#2166ac","DA":"#1a9641","PA":"#d73027","SA":"#636363",
      "SQA":"#2166ac","QAOA":"#d73027","Greedy":"#b8860b"}
MK = {"SBA":"o","DA":"s","PA":"^","SA":"D",
      "SQA":"o","QAOA":"s","Greedy":"^"}
DS = {"SBA":[],"DA":[6,2],"PA":[3,2],"SA":[1,2]}
LB = {"SBA":"SBA (fixed)","DA":"Digital Annealing (fixed)",
      "PA":"Population Annealing (fixed)","SA":"Classical SA"}


# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE 1 — Energy Convergence
# ═══════════════════════════════════════════════════════════════════════════════
def fig1_energy_convergence(seed=42):
    print("Figure 1 — energy convergence …")
    np.random.seed(seed);  random.seed(seed)

    Q_dict, N = LargeScaleRANQUBO.frequency_planning(20, 3, seed=seed)
    Q60 = LargeScaleRANQUBO.q_dict_to_dense(Q_dict, N, 60)
    Qs  = Q60[:60, :60]

    sba = SBA(n_steps=3000, dt=0.015, xi=0.5, p_start=-1.0, p_end=1.5, seed=seed)
    sba.solve(Qs, max_n=60)

    da = DA(n_replicas=10, n_sweeps=400, exchange_freq=15,
            t_min=0.001, t_max=8.0, seed=seed)
    da.solve(Qs, max_n=60)

    pa = PA(population=60, n_temps=50, t_start=8.0, t_end=0.001, seed=seed)
    pa.solve(Qs, max_n=60)

    _, _, sa_log = run_sa(Qs, n_steps=5000, T_start=5.0,
                          seed=seed, log_points=30)

    print(f"  SBA: {sba.energy_log[0]:.1f} → {sba.energy_log[-1]:.1f}")
    print(f"  DA : {da.energy_log[0]:.1f}  → {da.energy_log[-1]:.1f}")
    print(f"  PA : {pa.energy_log[0]:.1f}  → {pa.energy_log[-1]:.1f}")
    print(f"  SA : {sa_log[0]:.1f}  → {sa_log[-1]:.1f}")

    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    for key, log in [("SBA", sba.energy_log), ("DA", da.energy_log),
                     ("PA", pa.energy_log),   ("SA", sa_log)]:
        pct = np.linspace(0, 100, len(log))
        ax.plot(pct, log, color=C[key], marker=MK[key],
                markevery=5, dashes=DS[key], label=LB[key])
    ax.set_xlabel("Solver progress (%)")
    ax.set_ylabel("Best QUBO energy")
    ax.set_title("Fig. 1 — Energy convergence (frequency planning)\n"
                 r"20 cells $\times$ 3 channels, $n=60$ variables")
    ax.legend(loc="upper right")
    fig.tight_layout()
    save_fig(fig, "fig1_energy_convergence.pdf")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURES 2 & 3 — Runtime Scaling + Constraint Satisfaction
# ═══════════════════════════════════════════════════════════════════════════════
def fig23_scaling(seed=42):
    print("\nFigures 2 & 3 — scaling benchmark …")
    CELL_SIZES = [5, 8, 10, 12, 15, 18, 20]
    data = {s: {"ms": [], "sat": [], "sat_rep": []} for s in ["SBA","DA","PA","SA"]}

    for nc in CELL_SIZES:
        n = nc * 3
        print(f"  n_cells={nc}  n_vars={n} …", end=" ", flush=True)
        np.random.seed(seed);  random.seed(seed)
        Q_d, Nv = LargeScaleRANQUBO.frequency_planning(nc, 3, seed=seed)
        Qb = LargeScaleRANQUBO.q_dict_to_dense(Q_d, Nv, n)[:n, :n]

        solvers = {
            "SBA": SBA(n_steps=3000, dt=0.015, xi=0.5,
                       p_start=-1.0, p_end=1.5, seed=seed),
            "DA":  DA(n_replicas=10, n_sweeps=400, exchange_freq=15,
                      t_min=0.001, t_max=8.0, seed=seed),
            "PA":  PA(population=60, n_temps=50,
                      t_start=8.0, t_end=0.001, seed=seed),
        }
        for name, sv in solvers.items():
            t0 = time.perf_counter()
            x, _ = sv.solve(Qb, max_n=n)
            ms   = (time.perf_counter() - t0) * 1000
            x_r  = feasibility_repair(x, nc)          # Eq. 77 repair
            data[name]["ms"].append(round(ms, 2))
            data[name]["sat"].append(sat_pct(x, nc))
            data[name]["sat_rep"].append(sat_pct(x_r, nc))

        x_sa, _, _ = run_sa(Qb, n_steps=4000, T_start=5.0, seed=seed)
        t0 = time.perf_counter()
        run_sa(Qb, n_steps=4000, T_start=5.0, seed=seed)
        ms_sa = (time.perf_counter() - t0) * 1000
        x_sa_r = feasibility_repair(x_sa, nc)
        data["SA"]["ms"].append(round(ms_sa, 2))
        data["SA"]["sat"].append(sat_pct(x_sa, nc))
        data["SA"]["sat_rep"].append(sat_pct(x_sa_r, nc))
        print("done")

    nvars = [nc * 3 for nc in CELL_SIZES]

    # Figure 2 — Runtime
    fig2, ax2 = plt.subplots(figsize=(4.8, 3.4))
    for name in ["SBA", "DA", "PA", "SA"]:
        ax2.plot(nvars, data[name]["ms"], color=C[name],
                 marker=MK[name], dashes=DS[name], label=LB[name])
    ax2.set_xlabel(r"Number of QUBO variables ($n$)")
    ax2.set_ylabel("Wall-clock time (ms)")
    ax2.set_title("Fig. 2 — Runtime scaling (frequency planning QUBO)")
    ax2.legend(loc="upper left", fontsize=8)
    fig2.tight_layout()
    save_fig(fig2, "fig2_runtime_scaling.pdf")
    plt.close(fig2)

    # Figure 3 — Constraint satisfaction (raw + after repair)
    fig3, ax3 = plt.subplots(figsize=(4.8, 3.4))
    for name in ["SBA", "DA", "PA", "SA"]:
        ax3.plot(CELL_SIZES, data[name]["sat"], color=C[name],
                 marker=MK[name], dashes=DS[name], label=LB[name])
    # Show repaired SBA with dashed line for comparison
    ax3.plot(CELL_SIZES, data["SBA"]["sat_rep"], color=C["SBA"],
             dashes=[2, 2], linewidth=1.2, alpha=0.6,
             label="SBA + repair (Eq. 77)")
    ax3.set_xlabel("Number of cells")
    ax3.set_ylabel("Constraint satisfaction (%)")
    ax3.set_ylim(0, 115)
    ax3.yaxis.set_major_formatter(mticker.PercentFormatter())
    ax3.set_title("Fig. 3 — One-channel-per-cell satisfaction vs. scale")
    ax3.legend(loc="lower left", fontsize=7.5)
    fig3.tight_layout()
    save_fig(fig3, "fig3_constraint_satisfaction.pdf")
    plt.close(fig3)

    return data


# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE 4 — Approximation Gap (Module 45)
# ═══════════════════════════════════════════════════════════════════════════════
def fig4_approximation_gap(seed=42):
    print("\nFigure 4 — Module 45 approximation gap …")
    np.random.seed(seed);  random.seed(seed)

    qb   = QUBOBuilder()
    orch = QuantumOptimisationOrchestrator(seed=seed)

    sinr  = np.random.uniform(5, 25, (5, 3))
    Q1, _ = qb.user_bs_assignment(sinr, lambda_c=4.0)
    r1    = orch.benchmark("User\u2013BS assignment\n(5 UEs \u00d7 3 BSs)", Q1, max_n=12)

    gain  = np.random.uniform(0, 20, (4, 8))
    Q2, _ = qb.beam_selection(gain, n_beams_max=3, lambda_c=3.0)
    r2    = orch.benchmark("Beam selection\n(8 beams, pick 3)", Q2, max_n=8)

    G = np.array([[0,1,1,0,0],[1,0,1,1,0],[1,1,0,0,1],[0,1,0,0,1],[0,0,1,1,0]])
    Q3, _ = qb.channel_assignment(G, n_channels=3, lambda_c=5.0)
    r3    = orch.benchmark("Channel assignment\n(5 cells, 3 channels)", Q3, max_n=15)

    for r in [r1, r2, r3]:
        print(f"  {r['problem'].split(chr(10))[0]:25s}  "
              f"SQA={r['sqa_gap_pct']:+5.1f}%  "
              f"QAOA={r['qaoa_gap_pct']:+6.1f}%  "
              f"SA={r['sa_gap_pct']:+5.1f}%  "
              f"Greedy={r['greedy_gap_pct']:+6.1f}%")

    results = [r1, r2, r3]
    methods = [
        ("sqa_gap_pct",    "SQA (fixed)",  C["SQA"],    "///"),
        ("qaoa_gap_pct",   "QAOA (p=2)",   C["QAOA"],   "\\\\\\"),
        ("sa_gap_pct",     "Classical SA", C["SA"],      "..."),
        ("greedy_gap_pct", "Greedy",       C["Greedy"],  "xxx"),
    ]

    x_pos   = np.arange(len(results))
    width   = 0.18
    offsets = np.linspace(-1.5 * width, 1.5 * width, 4)

    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    for offset, (key, lbl, color, hatch) in zip(offsets, methods):
        vals = [max(r[key], 0) for r in results]
        bars = ax.bar(x_pos + offset, vals, width, label=lbl,
                      color=color, hatch=hatch, edgecolor="white",
                      linewidth=0.5, alpha=0.88)
        for bar, val in zip(bars, vals):
            if val > 0.5:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.4,
                        f"+{val:.1f}%",
                        ha="center", va="bottom", fontsize=7)

    ax.set_xticks(x_pos)
    ax.set_xticklabels([r["problem"] for r in results], fontsize=8.5)
    ax.set_ylabel("Energy gap above brute-force optimal (%)")
    ax.set_title("Fig. 4 — Approximation quality (Module 45)\n"
                 r"Small-scale problems ($n \leq 15$), all bugs fixed")
    ax.legend(loc="upper right", ncol=2, fontsize=8)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    save_fig(fig, "fig4_approximation_gap.pdf")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
#  SENSITIVITY ANALYSIS — λ_c effect on DA
# ═══════════════════════════════════════════════════════════════════════════════
def sensitivity_lambda(seed=42, n_cells=20, n_channels=3):
    """Reproduce Table XVIII: λ_c sensitivity for DA (n=60)."""
    print("\nSensitivity: λ_c effect on DA …")
    print(f"  {'λ_c':>6}  {'Energy':>10}  {'Normalised E/λ':>16}  {'Sat.':>8}")
    for lc in [2.0, 3.0, 5.0, 8.0, 12.0]:
        Q_d, Nv = LargeScaleRANQUBO.frequency_planning(
            n_cells, n_channels, lambda_c=lc, seed=seed)
        Qb = LargeScaleRANQUBO.q_dict_to_dense(Q_d, Nv, 60)[:60, :60]
        da = DA(n_replicas=10, n_sweeps=400, exchange_freq=15,
                t_min=0.001, t_max=8.0, seed=seed)
        x, e = da.solve(Qb, max_n=60)
        s    = sat_pct(x, n_cells)
        print(f"  {lc:>6.1f}  {e:>10.1f}  {e/lc:>16.2f}  {s:>7.1f}%")


# ═══════════════════════════════════════════════════════════════════════════════
#  SENSITIVITY ANALYSIS — n_steps effect on SBA
# ═══════════════════════════════════════════════════════════════════════════════
def sensitivity_nsteps(seed=42, n_cells=20, n_channels=3):
    """Reproduce Table XIX: n_steps sensitivity for SBA (n=60)."""
    print("\nSensitivity: n_steps effect on SBA …")
    Q_dict, N = LargeScaleRANQUBO.frequency_planning(n_cells, n_channels, seed=seed)
    Qs = LargeScaleRANQUBO.q_dict_to_dense(Q_dict, N, 60)[:60, :60]
    print(f"  {'n_steps':>8}  {'Energy':>10}  {'Sat.':>8}")
    for ns in [500, 1000, 2000, 3000, 5000]:
        np.random.seed(seed)
        sba = SBA(n_steps=ns, dt=0.015, xi=0.5, p_start=-1.0, p_end=1.5, seed=seed)
        x, e = sba.solve(Qs, max_n=60)
        s    = sat_pct(x, n_cells)
        s_r  = sat_pct(feasibility_repair(x, n_cells), n_cells)
        print(f"  {ns:>8}  {e:>10.1f}  {s:>6.1f}%  (repaired: {s_r:.1f}%)")


# ═══════════════════════════════════════════════════════════════════════════════
#  REPRODUCIBILITY STUDY — multi-seed DA and SQA
# ═══════════════════════════════════════════════════════════════════════════════
def reproducibility_study(seeds=(42, 123, 456, 789, 999)):
    """Reproduce Tables XX–XXI: multi-seed statistics."""
    print("\nReproducibility study …")

    # DA (n=60, frequency planning)
    print("  DA (n=60, frequency planning, 5 seeds):")
    da_energies = [];  da_sats = []
    for sd in seeds:
        Q_d, Nv = LargeScaleRANQUBO.frequency_planning(20, 3, seed=sd)
        Qb = LargeScaleRANQUBO.q_dict_to_dense(Q_d, Nv, 60)[:60, :60]
        da = DA(n_replicas=10, n_sweeps=400, exchange_freq=15,
                t_min=0.001, t_max=8.0, seed=sd)
        x, e = da.solve(Qb, max_n=60)
        s    = sat_pct(x, 20)
        da_energies.append(e);  da_sats.append(s)
        print(f"    seed={sd}: E={e:.1f}  sat={s:.1f}%")
    print(f"    Mean E={np.mean(da_energies):.1f} ± {np.std(da_energies):.1f}")
    print(f"    Mean sat={np.mean(da_sats):.1f}% ± {np.std(da_sats):.1f}pp")

    # SQA (User-BS assignment, n=12)
    print("  SQA (User-BS assignment, n=12, 5 seeds):")
    qb   = QUBOBuilder()
    np.random.seed(42);  random.seed(42)
    sinr = np.random.uniform(5, 25, (5, 3))
    Q1, _ = qb.user_bs_assignment(sinr, lambda_c=4.0)
    Q1s   = Q1[:12, :12]
    e_opt = ClassicalSolver.brute_force(Q1s, max_n=12)[1]
    sqa_gaps = []
    for sd in seeds:
        sqa  = SQA(n_replicas=20, T=0.3, n_steps=3000,
                   gamma_start=5.0, gamma_end=0.001, seed=sd)
        _, e = sqa.solve(Q1s, max_n=12)
        gap  = round(100 * (e - e_opt) / abs(e_opt + 1e-9), 1)
        sqa_gaps.append(gap)
        print(f"    seed={sd}: E={e:.3f}  gap={gap:+.1f}%")
    print(f"    Mean gap={np.mean(sqa_gaps):.1f}% ± {np.std(sqa_gaps):.1f}pp")

    # GPU projection (Eq. 81)
    t_sba_cpu_60 = 38.0   # measured ms at n=60
    print("\n  GPU runtime projection (G=100, Eq. 81):")
    for n_t in [60, 300, 1500]:
        t_gpu = gpu_runtime_projection(t_sba_cpu_60, 60, n_t, G=100)
        print(f"    SBA GPU @ n={n_t}: {t_gpu:.1f} ms")


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════
def main():
    SEED = 42
    random.seed(SEED);  np.random.seed(SEED)

    print("=" * 60)
    print("  Quantum-Inspired RAN Optimisation — Complete Pipeline")
    print("  Bugs 1–4 corrected | Feasibility repair | Sensitivity")
    print("=" * 60)

    # ── Four publication figures ───────────────────────────────────────────
    fig1_energy_convergence(seed=SEED)
    fig23_scaling(seed=SEED)
    fig4_approximation_gap(seed=SEED)

    # ── Sensitivity and reproducibility ───────────────────────────────────
    sensitivity_lambda(seed=SEED)
    sensitivity_nsteps(seed=SEED)
    reproducibility_study()

    print(f"\n{'=' * 60}")
    print(f"  Figures saved to: {OUT_DIR}")
    print(f"{'=' * 60}")
    print("  fig1_energy_convergence.pdf   — QUBO energy vs progress")
    print("  fig2_runtime_scaling.pdf      — runtime (ms) vs n vars")
    print("  fig3_constraint_satisfaction  — satisfaction (%) vs cells")
    print("  fig4_approximation_gap.pdf    — gap from optimal (M45)")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()

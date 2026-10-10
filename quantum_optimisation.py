"""
=============================================================================
  MODULE 45 — Quantum-Assisted Optimisation for RAN
=============================================================================
Many RAN optimisation problems are NP-hard: user scheduling, beam assignment,
slice admission, frequency planning. Quantum computing (and quantum-inspired
classical algorithms) can attack these more efficiently.

This module simulates quantum-assisted approaches using:
  • QUBO (Quadratic Unconstrained Binary Optimisation) formulation
  • Simulated Quantum Annealing (SQA) — classical simulation of QA
  • QAOA (Quantum Approximate Optimisation Algorithm) — variational circuit
  • D-Wave-style Ising model for combinatorial RAN problems
  • Benchmark: classical vs quantum-inspired vs brute-force

RAN problems formulated as QUBO:
  1. User-BS assignment (max-SINR with interference constraints)
  2. Beam selection (M beams, K UEs → K×M binary assignment)
  3. Network slicing admission (0/1 per slice request)
  4. Frequency channel assignment (graph colouring)
=============================================================================
"""

import math
import random
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Callable
from itertools import product


# ─────────────────────────────────────────────────────────────────────────────
# 45.1  QUBO formulation for RAN problems
# ─────────────────────────────────────────────────────────────────────────────

class QUBOBuilder:
    """
    Translates RAN combinatorial problems into QUBO form:
      minimise  x^T Q x     (x ∈ {0,1}^n)

    The matrix Q encodes:
      - Objective: reward for good assignments (negative diagonal)
      - Constraints: penalty for violations (large positive off-diagonal)
    """

    @staticmethod
    def user_bs_assignment(
        sinr_matrix:   np.ndarray,   # [K × M] SINR (dB) of UE k to BS m
        max_load:      float = 0.8,  # max fractional BS load
        lambda_c:      float = 5.0,  # constraint penalty weight
    ) -> Tuple[np.ndarray, int]:
        """
        Assign K UEs to M BSs to maximise total SINR.
        Binary variable x_{km} = 1 iff UE k assigned to BS m.

        Constraints:
          (1) Each UE assigned to exactly one BS: Σ_m x_{km} = 1  ∀k
          (2) BS load limit: Σ_k x_{km} ≤ max_load·K/M  ∀m

        Returns (Q matrix [KM × KM], n_vars).
        """
        K, M  = sinr_matrix.shape
        N     = K * M   # total binary variables
        Q     = np.zeros((N, N))

        def idx(k, m): return k * M + m

        # Objective: maximise SINR → minimise -SINR
        snr_lin = 10 ** (sinr_matrix / 10)
        snr_norm = snr_lin / (snr_lin.max() + 1e-9)
        for k in range(K):
            for m in range(M):
                Q[idx(k,m), idx(k,m)] -= snr_norm[k, m]

        # Constraint 1: each UE assigned to exactly one BS
        # ||Σ_m x_{km} - 1||² = Σ_m x_{km}² + 2Σ_{m<m'} x_{km}x_{km'} - 2Σ_m x_{km} + 1
        for k in range(K):
            for m in range(M):
                Q[idx(k,m), idx(k,m)] += lambda_c * (1 - 2)  # x² - 2x term
            for m in range(M):
                for mp in range(m+1, M):
                    Q[idx(k,m), idx(k,mp)] += 2 * lambda_c    # cross term

        # Constraint 2: load per BS (soft penalty for over-load)
        cap = max(1, int(max_load * K / M))
        for m in range(M):
            for k in range(K):
                for kp in range(k+1, K):
                    # Penalty if both k and k' go to m and load exceeded
                    if 2 > cap:
                        Q[idx(k,m), idx(kp,m)] += lambda_c * 0.5

        return Q, N

    @staticmethod
    def beam_selection(
        gain_matrix: np.ndarray,   # [K × B] beam gain for UE k, beam b
        n_beams_max: int = 4,      # max simultaneous beams
        lambda_c:    float = 3.0,
    ) -> Tuple[np.ndarray, int]:
        """
        Select up to n_beams_max beams from B candidates for K UEs.
        x_b ∈ {0,1}: whether beam b is activated.
        Maximise total gain subject to max beams constraint.
        """
        K, B = gain_matrix.shape
        Q    = np.zeros((B, B))

        # Objective: total gain from activated beams
        total_gain = gain_matrix.sum(axis=0)  # [B] gain summed over UEs
        gain_norm  = total_gain / (total_gain.max() + 1e-9)
        for b in range(B):
            Q[b, b] -= gain_norm[b]

        # Constraint: at most n_beams_max beams active
        # Penalty: (Σ_b x_b - n_max)² if Σ > n_max (quadratic slack)
        for b in range(B):
            for bp in range(b+1, B):
                Q[b, bp] += 2 * lambda_c / (n_beams_max * B)

        return Q, B

    @staticmethod
    def channel_assignment(
        interference_graph: np.ndarray,   # [N_cells × N_cells] adjacency
        n_channels:         int = 3,
        lambda_c:           float = 4.0,
    ) -> Tuple[np.ndarray, int]:
        """
        Graph colouring: assign one of n_channels to each cell
        such that adjacent cells have different channels.
        x_{ic} = 1 iff cell i uses channel c.
        """
        N_cells = len(interference_graph)
        N       = N_cells * n_channels
        Q       = np.zeros((N, N))

        def idx(i, c): return i * n_channels + c

        # Each cell gets exactly one channel
        for i in range(N_cells):
            for c in range(n_channels):
                Q[idx(i,c), idx(i,c)] += lambda_c * (1 - 2)
            for c in range(n_channels):
                for cp in range(c+1, n_channels):
                    Q[idx(i,c), idx(i,cp)] += 2 * lambda_c

        # Adjacent cells must use different channels
        for i in range(N_cells):
            for j in range(i+1, N_cells):
                if interference_graph[i, j] > 0:
                    for c in range(n_channels):
                        Q[idx(i,c), idx(j,c)] += lambda_c * 2

        return Q, N


# ─────────────────────────────────────────────────────────────────────────────
# 45.2  Simulated Quantum Annealing (SQA)
# ─────────────────────────────────────────────────────────────────────────────

class SimulatedQuantumAnnealer:
    """
    Simulated Quantum Annealing (SQA) — classical simulation of quantum annealing.

    Introduces a transverse-field term Γ that acts as quantum tunnelling:
    H(s) = -(1-s) Σ_i σ^x_i  +  s · H_problem(σ^z)
    where s ∈ [0,1] increases from 0 (full quantum) to 1 (classical).

    Implemented via path-integral Monte Carlo with P replicas (Trotter slices).
    This is not the same as running on actual quantum hardware but captures
    the tunnelling-assisted escape from local minima.
    """

    def __init__(
        self,
        n_replicas:  int   = 8,     # Trotter slices (more = better simulation)
        T:           float = 0.5,   # temperature for thermal fluctuations
        n_steps:     int   = 500,   # annealing steps
        gamma_start: float = 3.0,   # initial transverse field
        gamma_end:   float = 0.01,  # final transverse field
        seed:        int   = 42,
    ):
        self.P           = n_replicas
        self.T           = T
        self.n_steps     = n_steps
        self.gamma_start = gamma_start
        self.gamma_end   = gamma_end
        self.rng         = np.random.default_rng(seed)
        self.energy_log: List[float] = []

    def _qubo_to_ising(self, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Convert QUBO (x ∈ {0,1}) to Ising (σ ∈ {-1,+1}).
        x_i = (σ_i + 1)/2
        E_ising = Σ_{ij} J_{ij} σ_i σ_j + Σ_i h_i σ_i + const
        """
        N = Q.shape[0]
        J = np.zeros((N, N))
        h = np.zeros(N)
        for i in range(N):
            h[i] = Q[i,i]/2 + sum(Q[i,j]/4 + Q[j,i]/4 for j in range(N) if j != i)
            for j in range(i+1, N):
                J[i,j] = (Q[i,j] + Q[j,i]) / 4
        return J, h

    def _ising_energy(
        self,
        spins: np.ndarray,   # [N] ∈ {-1,+1}
        J:     np.ndarray,   # [N×N]
        h:     np.ndarray,   # [N]
    ) -> float:
        E = 0.0
        N = len(spins)
        for i in range(N):
            E += h[i] * spins[i]
            for j in range(i+1, N):
                E += J[i,j] * spins[i] * spins[j]
        return E

    def solve(
        self,
        Q:        np.ndarray,
        max_n:    int = 20,   # cap problem size for speed
    ) -> Tuple[np.ndarray, float]:
        """
        Solve QUBO using SQA. Returns (best_x, best_energy).
        """
        N   = min(Q.shape[0], max_n)
        Qs  = Q[:N, :N]
        J, h = self._qubo_to_ising(Qs)

        # Initialise P replicas randomly
        replicas = self.rng.choice([-1, 1], size=(self.P, N))

        best_spins  = replicas[0].copy()
        best_energy = self._ising_energy(replicas[0], J, h)

        for step in range(self.n_steps):
            # Anneal schedule
            s       = step / self.n_steps
            gamma   = self.gamma_start * (1-s) + self.gamma_end * s
            beta_eff= 1.0 / max(self.T, 1e-9)

            # Quantum coupling between replicas (transverse field)
            J_perp  = -0.5 * self.T * math.log(math.tanh(
                gamma * beta_eff / self.P)) if self.P > 1 else 0.0

            # Update each replica with Metropolis
            for p in range(self.P):
                i = self.rng.integers(0, N)
                # Local field from Ising couplings
                dE_local = 2 * replicas[p, i] * (
                    h[i] + sum(J[min(i,j), max(i,j)] * replicas[p, j]
                                for j in range(N) if j != i))
                # Quantum coupling to neighbouring replicas
                p_prev = (p - 1) % self.P
                p_next = (p + 1) % self.P
                dE_q   = 2 * J_perp * replicas[p, i] * (
                    replicas[p_prev, i] + replicas[p_next, i])
                dE     = dE_local + dE_q

                if dE < 0 or self.rng.random() < math.exp(-dE * beta_eff):
                    replicas[p, i] *= -1

                # Track best
                E = self._ising_energy(replicas[p], J, h)
                if E < best_energy:
                    best_energy = E
                    best_spins  = replicas[p].copy()

            if step % 50 == 0:
                self.energy_log.append(best_energy)

        # Convert back to QUBO x ∈ {0,1}
        best_x = np.zeros(Q.shape[0], dtype=int)
        best_x[:N] = ((best_spins + 1) // 2).astype(int)

        # QUBO energy
        qubo_energy = float(best_x @ Q @ best_x)
        return best_x, qubo_energy


# ─────────────────────────────────────────────────────────────────────────────
# 45.3  QAOA (variational quantum circuit simulation)
# ─────────────────────────────────────────────────────────────────────────────

class QAOASimulator:
    """
    Quantum Approximate Optimisation Algorithm (QAOA) simulation.

    QAOA uses p layers of alternating unitaries:
      U(β,γ) = e^{-iβ B} e^{-iγ C}
    where C is the problem Hamiltonian and B = Σ_i σ^x_i is the mixer.

    For small problems (n ≤ 16), we can exactly simulate the quantum state.
    For larger problems we use a classical approximation.

    The expectation value ⟨C⟩ is optimised over (β, γ) parameters.
    """

    def __init__(self, p_layers: int = 2, seed: int = 0):
        self.p   = p_layers
        self.rng = np.random.default_rng(seed)
        self.opt_log: List[float] = []

    def _problem_expectation_approx(
        self,
        Q:     np.ndarray,
        gamma: np.ndarray,
        beta:  np.ndarray,
    ) -> float:
        """
        Classical approximation of QAOA expectation value.
        For degree-2 QUBO, the p=1 analytical formula is:
        ⟨C⟩ ≈ Σ_{ij} Q_{ij}/2 · f(γ,β,graph_neighbourhood)
        We use a simplified bound.
        """
        N     = Q.shape[0]
        total = 0.0
        # Approximate: each variable contributes Q_{ii} modulated by γ,β
        for i in range(N):
            # Single-qubit contribution
            qi = Q[i,i]
            phase = sum(gamma) * qi
            total += qi * math.cos(phase) * math.sin(2 * sum(beta))
            # Two-qubit contributions
            for j in range(i+1, N):
                qij = (Q[i,j] + Q[j,i]) / 2
                phase2 = sum(gamma) * qij
                total += qij * math.sin(phase2) * math.cos(2 * sum(beta))
        return total / max(N, 1)

    def optimise_angles(
        self,
        Q:       np.ndarray,
        n_iter:  int = 40,
    ) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Optimise QAOA angles (β, γ) via gradient-free search.
        Returns (best_beta, best_gamma, best_expectation).
        """
        best_val   = float('inf')
        best_beta  = np.zeros(self.p)
        best_gamma = np.zeros(self.p)

        for _ in range(n_iter):
            gamma = self.rng.uniform(0, math.pi, self.p)
            beta  = self.rng.uniform(0, math.pi/2, self.p)
            val   = self._problem_expectation_approx(Q, gamma, beta)
            self.opt_log.append(val)
            if val < best_val:
                best_val   = val
                best_gamma = gamma.copy()
                best_beta  = beta.copy()

        return best_beta, best_gamma, best_val

    def sample_solution(
        self,
        Q:         np.ndarray,
        best_gamma:np.ndarray,
        best_beta: np.ndarray,
        n_shots:   int = 100,
    ) -> np.ndarray:
        """
        Sample candidate solutions using approximate QAOA distribution.
        Returns the best (lowest energy) binary string found.
        """
        N        = Q.shape[0]
        best_x   = None
        best_e   = float('inf')

        for _ in range(n_shots):
            # Approximate sampling: bias toward low-energy states
            probs = np.zeros(N)
            for i in range(N):
                theta_i = 2 * sum(best_beta) + Q[i,i] * sum(best_gamma)
                probs[i] = (1 + math.cos(theta_i)) / 2
            probs = np.clip(probs, 0.05, 0.95)
            x = (self.rng.random(N) < probs).astype(int)
            e = float(x @ Q @ x)
            if e < best_e:
                best_e = e
                best_x = x.copy()

        return best_x, best_e


# ─────────────────────────────────────────────────────────────────────────────
# 45.4  Classical baselines
# ─────────────────────────────────────────────────────────────────────────────

class ClassicalSolver:
    """
    Classical QUBO solvers for benchmarking.
    """

    @staticmethod
    def greedy(Q: np.ndarray) -> Tuple[np.ndarray, float]:
        """Greedy: set each variable to 1 if it reduces energy."""
        N   = Q.shape[0]
        x   = np.zeros(N, dtype=int)
        for i in range(N):
            e0 = float((x @ Q @ x))
            x[i] = 1
            e1 = float((x @ Q @ x))
            if e1 >= e0:
                x[i] = 0
        return x, float(x @ Q @ x)

    @staticmethod
    def simulated_annealing(
        Q:       np.ndarray,
        T_start: float = 2.0,
        T_end:   float = 0.01,
        n_steps: int   = 2000,
        seed:    int   = 0,
    ) -> Tuple[np.ndarray, float]:
        """Classical simulated annealing on QUBO."""
        rng  = np.random.default_rng(seed)
        N    = Q.shape[0]
        x    = rng.integers(0, 2, N)
        E    = float(x @ Q @ x)
        best_x, best_E = x.copy(), E

        for step in range(n_steps):
            T   = T_start * (T_end/T_start) ** (step/n_steps)
            i   = rng.integers(0, N)
            x[i] ^= 1
            E_new = float(x @ Q @ x)
            dE    = E_new - E
            if dE < 0 or rng.random() < math.exp(-dE/max(T,1e-9)):
                E = E_new
                if E < best_E:
                    best_E, best_x = E, x.copy()
            else:
                x[i] ^= 1

        return best_x, best_E

    @staticmethod
    def brute_force(Q: np.ndarray, max_n: int = 16) -> Tuple[np.ndarray, float]:
        """Exact brute-force (only feasible for n ≤ 20)."""
        N      = min(Q.shape[0], max_n)
        Qs     = Q[:N, :N]
        best_x = np.zeros(N, dtype=int)
        best_E = float('inf')
        for bits in range(2**N):
            x = np.array([(bits >> i) & 1 for i in range(N)], dtype=int)
            E = float(x @ Qs @ x)
            if E < best_E:
                best_E, best_x = E, x
        return best_x, best_E


# ─────────────────────────────────────────────────────────────────────────────
# 45.5  Quantum optimisation orchestrator
# ─────────────────────────────────────────────────────────────────────────────

class QuantumOptimisationOrchestrator:
    """
    Benchmarks quantum-inspired vs classical solvers on RAN QUBO problems.
    """

    def __init__(self, seed: int = 42):
        self.sqa    = SimulatedQuantumAnnealer(n_replicas=8, n_steps=500, seed=seed)
        self.qaoa   = QAOASimulator(p_layers=2, seed=seed)
        self.clas   = ClassicalSolver()
        self.results: List[Dict] = []

    def benchmark(
        self,
        problem_name: str,
        Q:            np.ndarray,
        max_n:        int = 18,
    ) -> Dict:
        """Run all solvers and compare."""
        N = min(Q.shape[0], max_n)
        Qs = Q[:N, :N]

        # SQA
        x_sqa, e_sqa   = self.sqa.solve(Qs, max_n=N)

        # QAOA
        bg, gm, _       = self.qaoa.optimise_angles(Qs, n_iter=30)
        x_qaoa, e_qaoa  = self.qaoa.sample_solution(Qs, gm, bg, n_shots=100)

        # Classical SA
        x_sa, e_sa      = self.clas.simulated_annealing(Qs, seed=42)

        # Greedy
        x_gr, e_gr      = self.clas.greedy(Qs)

        # Brute force (only if small)
        if N <= 16:
            x_bf, e_bf  = self.clas.brute_force(Qs, max_n=N)
        else:
            x_bf, e_bf  = None, None

        # Quality metrics (how close to optimal/best known)
        best_known = e_bf if e_bf is not None else min(e_sqa, e_sa, e_qaoa, e_gr)
        def gap(e):
            if best_known == 0:
                return 0.0
            return round(100 * (e - best_known) / abs(best_known + 1e-9), 1)

        result = {
            "problem":      problem_name,
            "n_vars":       N,
            "sqa_energy":   round(e_sqa, 3),
            "qaoa_energy":  round(e_qaoa, 3),
            "sa_energy":    round(e_sa, 3),
            "greedy_energy":round(e_gr, 3),
            "optimal_energy": round(e_bf, 3) if e_bf is not None else None,
            "sqa_gap_pct":  gap(e_sqa),
            "qaoa_gap_pct": gap(e_qaoa),
            "sa_gap_pct":   gap(e_sa),
            "greedy_gap_pct": gap(e_gr),
        }
        self.results.append(result)
        return result

    def print_result(self, r: Dict):
        print(f"\n  Problem: {r['problem']}  (n={r['n_vars']} variables)")
        print(f"  {'Method':<22}  {'Energy':>10}  {'Gap to opt':>11}")
        print("  " + "─"*46)
        for name, e_key, g_key in [
            ("Brute Force (optimal)", "optimal_energy", None),
            ("Simul. Quantum Anneal.", "sqa_energy",    "sqa_gap_pct"),
            ("QAOA (2 layers)",       "qaoa_energy",    "qaoa_gap_pct"),
            ("Classical SA",          "sa_energy",      "sa_gap_pct"),
            ("Greedy",               "greedy_energy",  "greedy_gap_pct"),
        ]:
            e = r.get(e_key)
            g = r.get(g_key)
            if e is None:
                continue
            gap_str = f"{g:+.1f}%" if g is not None else "  (optimal)"
            bar     = "▓" * max(0, 10 - int(abs(g or 0) / 2)) if g else "▓▓▓▓▓▓▓▓▓▓"
            print(f"  {name:<22}  {e:>10.3f}  {gap_str:>11}  {bar}")

    def summary(self):
        print("\n" + "="*65)
        print("  QUANTUM OPTIMISATION BENCHMARK — SUMMARY")
        print("="*65)
        print(f"  Problems solved: {len(self.results)}")
        for r in self.results:
            best = r["optimal_energy"]
            sqa_best = r["sqa_gap_pct"]
            print(f"\n  {r['problem']} (n={r['n_vars']})")
            print(f"    SQA gap:  {r['sqa_gap_pct']:+.1f}%   "
                  f"QAOA gap: {r['qaoa_gap_pct']:+.1f}%   "
                  f"SA gap: {r['sa_gap_pct']:+.1f}%   "
                  f"Greedy: {r['greedy_gap_pct']:+.1f}%")
        print("="*65)


# ─────────────────────────────────────────────────────────────────────────────
# Standalone demo
# ─────────────────────────────────────────────────────────────────────────────

def _demo():
    random.seed(42); np.random.seed(42)
    print("="*60)
    print("  MODULE 45 — Quantum-Assisted RAN Optimisation")
    print("="*60)

    orch = QuantumOptimisationOrchestrator(seed=42)
    qb   = QUBOBuilder()

    # Problem 1: User-BS assignment (6 UEs, 3 BSs)
    K, M = 6, 3
    sinr  = np.random.uniform(5, 25, (K, M))   # dB
    Q1, _ = qb.user_bs_assignment(sinr, max_load=0.8, lambda_c=4.0)
    r1    = orch.benchmark("User-BS Assignment (6×3)", Q1, max_n=12)
    orch.print_result(r1)

    # Problem 2: Beam selection (8 candidate beams, pick 3)
    gain  = np.random.uniform(0, 20, (4, 8))    # gain matrix [UEs × beams]
    Q2, _ = qb.beam_selection(gain, n_beams_max=3, lambda_c=3.0)
    r2    = orch.benchmark("Beam Selection (8 beams)", Q2, max_n=8)
    orch.print_result(r2)

    # Problem 3: Channel assignment / frequency planning (5 cells, 3 channels)
    G = np.array([                              # interference adjacency
        [0,1,1,0,0],
        [1,0,1,1,0],
        [1,1,0,0,1],
        [0,1,0,0,1],
        [0,0,1,1,0],
    ])
    Q3, _ = qb.channel_assignment(G, n_channels=3, lambda_c=5.0)
    r3    = orch.benchmark("Freq Channel Assignment", Q3, max_n=15)
    orch.print_result(r3)

    # Show QUBO structure insight
    print(f"\n  QUBO Q-matrix insights:")
    print(f"    User-BS ({K}×{M}): {Q1.shape[0]} binary vars, "
          f"density={np.count_nonzero(Q1)/Q1.size:.1%}")
    print(f"    Beam sel (8):    {Q2.shape[0]} binary vars, "
          f"density={np.count_nonzero(Q2)/Q2.size:.1%}")
    print(f"    Chan. assign:    {Q3.shape[0]} binary vars, "
          f"density={np.count_nonzero(Q3)/Q3.size:.1%}")

    # Annealing energy convergence
    print(f"\n  SQA energy convergence (User-BS problem):")
    if orch.sqa.energy_log:
        log = orch.sqa.energy_log
        for i, e in enumerate(log):
            bar = "█" * max(0, int(10 - (e - min(log)) / (max(log)-min(log)+1e-9) * 10))
            print(f"    step {i*50:>4}: E={e:>8.3f}  {bar}")

    orch.summary()


if __name__ == "__main__":
    _demo()

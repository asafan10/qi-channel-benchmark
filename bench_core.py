"""
bench_core.py -- additional checks (no new instance types).

Part A  (tuning parity)  native SA starting temperature T0 chosen on the held-out seeds 100-101 with the SAME grid
         and the same budget levels that were used to choose the polish temperature: T0 in {2, 0.6, 0.3, 0.15, 0.08},
         sweeps in {1, 4, 16, 64} x 1000, 3 runs each, K in {3, 4}, n = 200.  -> results/core_sa_t0.jsonl
Part B  (polish sensitivity) on the evaluation instances (seeds 0-9, n = 200, K in {3, 4}): one bSB batch per run at
         baseline steps, best-energy replica decoded, then native SA polish with 4000 sweeps for every T0 in the grid,
         4 runs.  -> results/core_pol_t0.jsonl
Part C  (independent-implementation diagnosis), seeds 0-4, n = 200, K in {3, 4}, 3 runs each:
         ours          : our bSB as in the paper (dt 0.5, 3000 steps, 8 replicas)
         ours-pkgset   : our bSB with the package's settings (dt 0.1, 3333 steps, 128 replicas, same linear pump)
         ours-ancilla  : our bSB with the linear term embedded as couplings to an extra spin, as the package does
                         (c0 then computed from the augmented matrix), paper settings
         pkg-ourset    : the package with our settings (time step 0.5, pressure slope 1/(0.5*3000), 3000 steps, 8 agents)
         -> results/core_sbdiag.jsonl
"""
import sys, time, math, json, warnings, numpy as np, scipy.sparse as sp
warnings.filterwarnings("ignore")
from exp_common import *

GRID = [2.0, 0.6, 0.3, 0.15, 0.08]


def part_a(out):
    done = done_keys(out, ("n", "seed", "K"))
    for K in (3, 4):
        for seed in (100, 101):
            if (200, seed, K) in done: continue
            W, csr, a, J, lam, jj = setup(200, seed, K); rng = np.random.default_rng(7000 + seed)
            res = {}
            for T0 in GRID:
                for f in (1, 4, 16, 64):
                    recs = []
                    for k in range(3):
                        t0 = time.perf_counter(); v = sa_native(rng.integers(0, K, 200), K, *csr, T0, 0.02, f * 1000, 50000 + 97 * k + f)[1]
                        recs.append(dict(f=float(v), t=time.perf_counter() - t0))
                    res[f"{T0}|{f}"] = recs
            append_jsonl(out, dict(n=200, seed=seed, K=K, runs=res)); print("A", K, seed, flush=True)


def part_b(out):
    done = done_keys(out, ("n", "seed", "K"))
    for K in (3, 4):
        for seed in range(10):
            if (200, seed, K) in done: continue
            W, csr, a, J, lam, jj = setup(200, seed, K); rng = np.random.default_rng(8100 + seed)
            res = {str(T0): [] for T0 in GRID}; fsel = []
            for b in range(4):
                Xb = sb_batch(a, J, CFG["sb_steps"], CFG["sb_R"], rng, False)
                es = [qubo_energy(x, a, *jj) for x in Xb]; sel = int(np.argmin(es))
                ch, _ = repair(Xb[sel].astype(np.int64), K, *csr); fsel.append(float(interference(ch, *csr)))
                for T0 in GRID:
                    t0 = time.perf_counter(); v = sa_native(ch, K, *csr, T0, 0.02, 4000, 61000 + 100 * b)[1]
                    res[str(T0)].append(dict(f=float(v), t=time.perf_counter() - t0))
            append_jsonl(out, dict(n=200, seed=seed, K=K, f_sel=fsel, runs=res)); print("B", K, seed, flush=True)


def ising_of(a, J):
    Jm = (J * 0.25).tocsr(); h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
    return Jm, h


def sb_ising(Jm, h, steps, R, rng, dt, ancilla=False):
    """ballistic SB on E = s.h + 0.5 s'Jm s ; with ancilla=True the field is embedded as couplings to spin N."""
    N = Jm.shape[0]
    if ancilla:
        hc = sp.csr_matrix(h.reshape(-1, 1))
        Jm = sp.bmat([[Jm, hc], [hc.T, None]]).tocsr(); h = np.zeros(N + 1); NN = N + 1
    else:
        NN = N
    c0 = 0.5 * math.sqrt(NN - 1) / (math.sqrt(Jm.multiply(Jm).sum()) + 1e-12)
    X = rng.uniform(-0.1, 0.1, (NN, R)); Y = rng.uniform(-0.1, 0.1, (NN, R))
    for t in range(steps):
        at = t / steps
        Y += dt * (-(1.0 - at) * X - c0 * (Jm @ X + h[:, None]))
        X += dt * Y
        m = np.abs(X) > 1; X[m] = np.sign(X[m]); Y[m] = 0.0
    S = np.where(X >= 0, 1, -1)
    if ancilla: S = S[:N] * S[N][None, :]
    return S.T                                          # R x N spins


def ising_energy(S, Jm, h):
    return S @ h + 0.5 * np.einsum("ri,ri->r", S, (Jm @ S.T).T)


def part_c(out):
    import torch, simulated_bifurcation as sbp
    torch.set_num_threads(1)
    done = done_keys(out, ("n", "seed", "K"))
    VARS = {"ours": dict(dt=0.5, steps=3000, R=8, anc=False), "ours-pkgset": dict(dt=0.1, steps=3333, R=128, anc=False),
            "ours-ancilla": dict(dt=0.5, steps=3000, R=8, anc=True)}
    for K in (3, 4):
        for seed in range(5):
            if (200, seed, K) in done: continue
            W, csr, a, J, lam, jj = setup(200, seed, K); Jm, h = ising_of(a, J); rng = np.random.default_rng(9200 + seed)
            res = {}
            for name, c in VARS.items():
                recs = []
                for k in range(3):
                    t0 = time.perf_counter(); S = sb_ising(Jm, h, c["steps"], c["R"], rng, c["dt"], c["anc"])
                    E = ising_energy(S, Jm, h); s = S[int(np.argmin(E))]; x = ((s + 1) // 2).astype(np.int64)
                    f, phi = decode_eval(x, K, csr); recs.append(dict(f=float(f), phi=float(phi), t=time.perf_counter() - t0))
                res[name] = recs
            Q = torch.tensor(Jm.toarray() / 2, dtype=torch.float32); l = torch.tensor(h, dtype=torch.float32)
            sbp.set_env(time_step=0.5, pressure_slope=1.0 / (0.5 * 3000)); recs = []
            for k in range(3):
                torch.manual_seed(5000 * seed + k); t0 = time.perf_counter()
                S, E = sbp.minimize(Q, l, domain="spin", agents=8, max_steps=3000, mode="ballistic", heated=False,
                                    verbose=False, best_only=True, early_stopping=False)
                x = ((S.cpu().numpy().astype(int) + 1) // 2).astype(np.int64); f, phi = decode_eval(x, K, csr)
                recs.append(dict(f=float(f), phi=float(phi), t=time.perf_counter() - t0))
            sbp.reset_env(); res["pkg-ourset"] = recs
            append_jsonl(out, dict(n=200, seed=seed, K=K, runs=res))
            print("C", K, seed, {k: (round(min(r["f"] for r in v), 1), round(np.mean([r["phi"] for r in v]), 2)) for k, v in res.items()}, flush=True)


if __name__ == "__main__":
    warm(3); warm(4)
    part_c("results/core_sbdiag.jsonl")
    part_a("results/core_sa_t0.jsonl")
    part_b("results/core_pol_t0.jsonl")
    print("DONE", flush=True)

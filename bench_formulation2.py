"""Formulation controls, part 2: separate the bit encoding from the move set.
qubo-swap : SA on the one-hot QUBO whose move switches the active bit of a cell to another channel (two bit flips evaluated
            together on the QUBO energy, so every state stays one-hot); rho in {0.25, 0.12, 0.06}
qubo-sameT: single-bit-flip SA on the QUBO, as in bench_formulation.py, adding rho = 0.06
Same schedule T 2.0 -> 0.02, same seeds 0-9, 200 cells, K = 3, two runs per setting. Output results/formulation2.jsonl."""
import math, time
import numpy as np
from numba import njit
from exp_common import *
from refs import build_ref
from bench_formulation import onehot


@njit(cache=True)
def sa_qubo_swap(x0, K, a, jp, jx, jd, T0, T1, sweeps, seed):
    np.random.seed(seed); N = len(x0); n = N // K; x = x0.copy(); g = a.copy()
    for i in range(N):
        for k in range(jp[i], jp[i + 1]):
            if x[jx[k]]: g[i] += jd[k]
    E = 0.0
    for i in range(N):
        if x[i]:
            E += a[i]
            for k in range(jp[i], jp[i + 1]):
                if jx[k] > i and x[jx[k]]: E += jd[k]
    best = E; bx = x.copy()
    for s in range(sweeps):
        T = T0 * (T1 / T0) ** (s / max(sweeps - 1, 1))
        for _ in range(n):
            cell = np.random.randint(n); c = -1
            for q in range(K):
                if x[cell * K + q]: c = q
            c2 = np.random.randint(K - 1)
            if c2 >= c: c2 += 1
            p = cell * K + c; q2 = cell * K + c2
            d1 = -g[p]                                   # switch p off
            jpq = 0.0
            for k in range(jp[q2], jp[q2 + 1]):
                if jx[k] == p: jpq = jd[k]
            d2 = g[q2] - jpq                             # then q2 on (field of q2 without p)
            d = d1 + d2
            if d <= 0 or np.random.random() < math.exp(-d / T):
                x[p] = 0; x[q2] = 1; E += d
                for k in range(jp[p], jp[p + 1]): g[jx[k]] -= jd[k]
                for k in range(jp[q2], jp[q2 + 1]): g[jx[k]] += jd[k]
        if E < best - 1e-12: best = E; bx = x.copy()
    return bx, best


if __name__ == "__main__":
    K, n = 3, 200; REF = build_ref(); out = "results/formulation2.jsonl"; warm(K)
    for seed in range(10):
        W, csr, *_ = setup(n, seed, K); rec = dict(n=n, seed=seed, K=K, ref=REF[(K, n, seed)], res={})
        for S in (4000, 16000, 64000):
            for r in range(2):
                sd = 1000 * seed + 10 * r + S % 7
                for rho in (0.25, 0.12, 0.06):
                    a2, J2, lam2 = build_qubo(W, K, rho); j2 = (J2.indptr.astype(np.int64), J2.indices.astype(np.int64), J2.data.astype(np.float64))
                    x0 = onehot(np.random.default_rng(sd).integers(0, K, n), K)
                    t0 = time.perf_counter(); bx, _ = sa_qubo_swap(x0, K, a2, *j2, 2.0, 0.02, S, sd)
                    rec["res"].setdefault(f"qubo-swap-rho{rho}|{S}", []).append(dict(f=float(decode_eval(bx, K, csr)[0]), t=time.perf_counter() - t0))
                    if rho == 0.06:
                        t0 = time.perf_counter(); bx, _ = sa_qubo(x0, a2, *j2, 2.0, 0.02, S, sd)
                        rec["res"].setdefault(f"qubo-sameT-rho0.06|{S}", []).append(dict(f=float(decode_eval(bx, K, csr)[0]), t=time.perf_counter() - t0))
        append_jsonl(out, rec); print("F2", seed, flush=True)

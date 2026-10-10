"""Formulation controls: is native SA's lead due to the algorithm or to the formulation?
Same annealer family, same geometric schedule T 2.0 -> 0.02 (objective units), same sweeps, K = 3, 200 cells, seeds 0-9:
  native      : assignment space, random cell recoloured (one-hot preserved by construction)       [native SA of the paper]
  native+empty: assignment space with an extra 'empty' state per cell costing lambda (the one-hot penalty of an empty cell)
  qubo-sameT  : single-bit-flip SA on the penalty QUBO with the native schedule (rho = 0.25 and 0.12); a QUBO sweep has nK proposals
  qubo-orig   : single-bit-flip SA on the QUBO with its own frozen schedule 0.5 lambda -> 0.15 (rho = 0.25), as in the paper
All final states decoded with the paper's repair map. Output results/formulation.jsonl."""
import json, math, sys, time
import numpy as np
from numba import njit
from exp_common import *
from refs import build_ref


@njit(cache=True)
def sa_empty(ch0, K, ip, ix, dt_, lam, T0, T1, sweeps, seed):
    """native SA with state K = empty (no interference, penalty lam)."""
    np.random.seed(seed); n = len(ch0); ch = ch0.copy(); E = 0.0
    for i in range(n):
        if ch[i] == K: E += lam
        for k in range(ip[i], ip[i + 1]):
            j = ix[k]
            if j > i and ch[i] == ch[j] and ch[i] < K: E += dt_[k]
    best = E; bch = ch.copy()
    for s in range(sweeps):
        T = T0 * (T1 / T0) ** (s / max(sweeps - 1, 1))
        for _ in range(n):
            i = np.random.randint(n); c2 = np.random.randint(K)
            if c2 >= ch[i]: c2 += 1
            d = 0.0
            if ch[i] == K: d -= lam
            if c2 == K: d += lam
            for k in range(ip[i], ip[i + 1]):
                cj = ch[ix[k]]
                if cj == K: continue
                if cj == c2: d += dt_[k]
                elif cj == ch[i]: d -= dt_[k]
            if d <= 0 or np.random.random() < math.exp(-d / T):
                ch[i] = c2; E += d
        if E < best - 1e-12: best = E; bch = ch.copy()
    return bch, best


def onehot(ch, K):
    x = np.zeros(len(ch) * K, dtype=np.int64)
    for i, c in enumerate(ch):
        if c < K: x[i * K + c] = 1
    return x


if __name__ == "__main__":
    K, n = 3, 200; REF = build_ref(); SWEEPS = (4000, 16000, 64000); out = "results/formulation.jsonl"; warm(K)
    for seed in range(10):
        W, csr, a, J, lam, jj = setup(n, seed, K); rec = dict(n=n, seed=seed, K=K, ref=REF[(K, n, seed)], res={})
        def add(name, f, t): rec["res"].setdefault(name, []).append(dict(f=float(f), t=t))
        for S in SWEEPS:
            for r in range(2):
                sd = 1000 * seed + 10 * r + S % 7
                t0 = time.perf_counter(); v = sa_native(np.random.default_rng(sd).integers(0, K, n), K, *csr, 2.0, 0.02, S, sd)[1]; add(f"native|{S}", v, time.perf_counter() - t0)
                t0 = time.perf_counter(); bch, _ = sa_empty(np.random.default_rng(sd).integers(0, K, n), K, *csr, lam, 2.0, 0.02, S, sd)
                add(f"native+empty|{S}", decode_eval(onehot(bch, K), K, csr)[0], time.perf_counter() - t0)
                for rho in (0.25, 0.12):
                    a2, J2, lam2 = build_qubo(W, K, rho); j2 = (J2.indptr.astype(np.int64), J2.indices.astype(np.int64), J2.data.astype(np.float64))
                    x0 = onehot(np.random.default_rng(sd).integers(0, K, n), K)
                    t0 = time.perf_counter(); bx, _ = sa_qubo(x0, a2, *j2, 2.0, 0.02, S, sd); add(f"qubo-sameT-rho{rho}|{S}", decode_eval(bx, K, csr)[0], time.perf_counter() - t0)
                x0 = onehot(np.random.default_rng(sd).integers(0, K, n), K)
                t0 = time.perf_counter(); bx, _ = sa_qubo(x0, a, *jj, 0.5 * lam, 0.15, S, sd); add(f"qubo-orig|{S}", decode_eval(bx, K, csr)[0], time.perf_counter() - t0)
        append_jsonl(out, rec); print("F", seed, flush=True)

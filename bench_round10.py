"""The constraint-preserving bSB of bench_round9.py compiled with Numba (same dynamics, same tuned
step size and coupling scale from round9_CPT), so that matched wall-clock comparisons are not handicapped by the NumPy implementation.
 NB    check against the NumPy version (K = 3 and 4, seeds 0-4, 3,000 and 12,000 steps, best plan of a batch of eight) and time per run.
 PL    planted 3-colourings (the instances of bench_round8.py part P and bench_round8q.py: every edge between cells of equal planted colour
       removed, optimum 0) at 50, 100 and 200 cells, five instances each: 32 replicas (four batches of eight) at 3,000 steps, as for bSB in
       Table S12, and at 200 cells also 12,000 and 48,000 steps.
 CPPN  matched wall-clock run pools (K = 3 or 4, seeds a-b), as CPP of round 9: 3,000, 12,000 and 48,000 steps, pools sized for 30, 100, 100 s.
Usage: python3 bench_round10.py NB|PL|CPPN a-b K  -> results/round10_<part>.json(l)"""
import sys, os, json, math, time, numpy as np
from numba import njit
from exp_common import *

part = sys.argv[1]; out = {}
CFG9 = {int(k): v for k, v in json.load(open("results/round9_CPP_cfg.json"))["cfg"].items()}


@njit(cache=True)
def _proj(v, K, u):
    for k in range(K): u[k] = v[k]
    for k in range(1, K):                       # insertion sort, descending (K <= 4)
        x = u[k]; j = k - 1
        while j >= 0 and u[j] < x:
            u[j + 1] = u[j]; j -= 1
        u[j + 1] = x
    css = 0.0; theta = 0.0
    for k in range(K):
        css += u[k]; t = (css - 1.0) / (k + 1)
        if u[k] - t > 0: theta = t
    for k in range(K):
        v[k] = max(v[k] - theta, 0.0)


@njit(cache=True)
def cp_bsb_nb(ip, ix, wd, n, K, steps, R, dt, c0, seed):
    np.random.seed(seed)
    out = np.zeros((R, n), np.int64); F = np.zeros((n, K)); v = np.zeros(K); u = np.zeros(K)
    for r in range(R):
        X = np.zeros((n, K)); Y = np.zeros((n, K))
        for i in range(n):
            for c in range(K): v[c] = 1.0 / K + np.random.uniform(-0.05, 0.05)
            _proj(v, K, u)
            m = 0.0
            for c in range(K): X[i, c] = v[c]; Y[i, c] = np.random.uniform(-0.05, 0.05); m += Y[i, c]
            for c in range(K): Y[i, c] -= m / K
        for t in range(steps):
            a = t / steps
            for i in range(n):
                for c in range(K): F[i, c] = -(1.0 - a) * (X[i, c] - 1.0 / K)
                for k in range(ip[i], ip[i + 1]):
                    j = ix[k]; w = wd[k]
                    for c in range(K): F[i, c] -= c0 * w * X[j, c]
            for i in range(n):
                m = 0.0
                for c in range(K): m += F[i, c]
                for c in range(K):
                    Y[i, c] += dt * (F[i, c] - m / K); v[c] = X[i, c] + dt * Y[i, c]
                _proj(v, K, u)
                for c in range(K):
                    X[i, c] = v[c]
                    if v[c] <= 0.0: Y[i, c] = 0.0
        for i in range(n):
            b = 0
            for c in range(1, K):
                if X[i, c] > X[i, b]: b = c
            out[r, i] = b
    return out


def run_batch(W, csr, K, steps, R, seed):
    n = W.shape[0]; N = n * K; cfg = CFG9[K]
    c0 = cfg["xi"] * 0.5 * math.sqrt(N - 1) / (math.sqrt(K) * math.sqrt((W ** 2).sum()))      # ||W kron I_K||_F = sqrt(K) ||W||_F
    C = cp_bsb_nb(csr[0], csr[1], csr[2], n, K, steps, R, cfg["dt"], c0, seed)
    return [float(interference(c, *csr)) for c in C]


def warm_cp():
    W = make_instance(12, 0, 3); csr = csr_from_W(W); run_batch(W, csr, 3, 5, 2, 0)


if part == "NB":
    from refs import build_ref
    REF = build_ref(); warm(3); warm(4); warm_cp(); rows = []
    for K in (3, 4):
        for seed in range(5):
            W = make_instance(200, seed, K); csr = csr_from_W(W); rec = dict(K=K, seed=seed, runs={})
            for steps in (3000, 12000):
                fs = []; ts = []
                for b in range(4):
                    t0 = time.perf_counter(); fs.append(min(run_batch(W, csr, K, steps, 8, 17000 + 100 * seed + 10 * b + K))); ts.append(time.perf_counter() - t0)
                rec["runs"][str(steps)] = dict(f=fs, t=float(np.mean(ts)))
            rec["ref"] = REF[(K, 200, seed)]; rows.append(rec)
            print("NB", K, seed, {s: (round((np.median(v["f"]) / rec["ref"] - 1) * 100, 2), round(v["t"], 3)) for s, v in rec["runs"].items()}, flush=True)
    out["rows"] = rows

if part == "PL":
    warm(3); warm_cp(); K = 3; rows = []
    for n in (50, 100, 200):
        for seed in range(5):
            W0 = make_instance(n, seed, K); col = np.random.default_rng(900 + seed).integers(0, K, n)
            W = W0 * (col[:, None] != col[None, :]); csr = csr_from_W(W)
            assert interference(col.astype(np.int64), *csr) == 0.0
            rec = dict(n=n, seed=seed, scale=float(W0[np.triu_indices(n, 1)].sum() / 3), res={})
            for steps in ((3000, 12000, 48000) if n == 200 else (3000,)):
                rec["res"][str(steps)] = [f for b in range(4) for f in run_batch(W, csr, K, steps, 8, 18000 + 100 * seed + 10 * b + n)]
            rows.append(rec); print("PL", n, seed, {k: (round(min(v), 3), int(sum(x == 0.0 for x in v))) for k, v in rec["res"].items()}, flush=True)
    out["rows"] = rows

if part == "CPPN":
    a_, b_ = map(int, sys.argv[2].split("-")); K = int(sys.argv[3]); warm(K); warm_cp(); jl = f"results/round10_CPPN_K{K}.jsonl"
    done_ = {json.loads(l)["seed"] for l in open(jl)} if os.path.exists(jl) else set()
    for seed in range(a_, b_):
        if seed in done_: continue
        W = make_instance(200, seed, K); csr = csr_from_W(W); rec = dict(K=K, seed=seed, cfg=CFG9[K], runs={}); k_ = 0
        for steps, Tmax in ((3000, 30), (12000, 100), (48000, 100)):
            pool = []; need = None
            while need is None or len(pool) < need:
                t0 = time.perf_counter(); f = min(run_batch(W, csr, K, steps, 8, 19000 + 1000 * seed + k_)); k_ += 1
                pool.append(dict(f=f, t=time.perf_counter() - t0))
                if need is None: need = max(4, min(150, 2 * int(Tmax // pool[0]["t"])))
            rec["runs"][str(steps)] = pool
        append_jsonl(jl, rec); print("CPPN", K, seed, {k: len(v) for k, v in rec["runs"].items()}, flush=True)
    sys.exit(0)

json.dump(out, open(f"results/round10_{part}.json", "w"), indent=1)

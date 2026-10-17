"""Single-bit-flip SA on the one-hot penalty QUBO (sa_qubo of bench_ran_qubo.py), tuned on held-out instances and
compared at matched wall-clock time on the instances of Table 7 (200 cells, K = 3 and 4, seeds 0-4), exactly as tuned SQA in bench_sqa_tuned.py.
 TUNE  grid on held-out seeds 100 and 101 at each K: penalty scale rho in {0.06, 0.12, 0.25}, starting temperature T0 in {0.15, 0.5, 1, 2,
       lambda/2} and final temperature T1 in {0.02, 0.005} (geometric schedule); 64,000 sweeps, two runs per seed, repaired cost scored against
       the better of two native SA runs of 64,000 sweeps (as for SQA).  -> results/saq_tune.json
 POOL  run pools with the tuned setting at 4,000, 16,000, 64,000 and 256,000 sweeps (twice the number of runs that fit in 30, 100, 100 and
       100 s, at most 150), two-process load.  -> results/saq_pool_K<K>.jsonl
Usage: python3 bench_saq_tuned.py TUNE | POOL a-b K"""
import sys, os, json, time, itertools, numpy as np
from exp_common import *

part = sys.argv[1]
GRID = list(itertools.product((0.06, 0.12, 0.25), (0.15, 0.5, 1.0, 2.0, "lam/2"), (0.02, 0.005)))


def saq_run(W, csr, K, cfg, sweeps, seed):
    a, J, lam = build_qubo(W, K, cfg["rho"]); jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64))
    T0 = 0.5 * lam if cfg["T0"] == "lam/2" else float(cfg["T0"])
    x0 = np.zeros(len(a), dtype=np.int64); ch = np.random.default_rng(seed).integers(0, K, W.shape[0]); x0[np.arange(W.shape[0]) * K + ch] = 1
    bx, _ = sa_qubo(x0, a, *jj, T0, cfg["T1"], sweeps, seed)
    f, phi = decode_eval(bx, K, csr)[:2]
    return float(f), float(phi)


if part == "TUNE":
    warm(3); warm(4); out = {}
    W0 = make_instance(20, 0, 3); a0, J0, _ = build_qubo(W0, 3, 0.25)
    sa_qubo(np.zeros(len(a0), dtype=np.int64), a0, J0.indptr.astype(np.int64), J0.indices.astype(np.int64), J0.data.astype(np.float64), 1.0, 0.1, 2, 1)
    for K in (3, 4):
        ref = {}
        for s in (100, 101):
            W = make_instance(200, s, K); csr = csr_from_W(W)
            ref[s] = min(sa_native(np.random.default_rng(s + k).integers(0, K, 200), K, *csr, 2.0, 0.02, 64000, 900 + k)[1] for k in range(2))
        res = {}
        for rho, T0, T1 in GRID:
            cfg = dict(rho=rho, T0=T0, T1=T1); v = []; ph = []; ts = []
            for s in (100, 101):
                W = make_instance(200, s, K); csr = csr_from_W(W)
                for k in range(2):
                    t0 = time.perf_counter(); f, phi = saq_run(W, csr, K, cfg, 64000, 33000 + 10 * s + k); ts.append(time.perf_counter() - t0)
                    v.append((f / ref[s] - 1) * 100); ph.append(phi)
            key = f"K{K}_rho{rho}_T0{T0}_T1{T1}"; res[key] = dict(cfg=cfg, mean_gap=float(np.mean(v)), phi=float(np.mean(ph)), t=float(np.mean(ts)))
            print("TUNE", key, round(float(np.mean(v)), 2), round(float(np.mean(ph)), 3), round(float(np.mean(ts)), 2), flush=True)
        best = min(res.values(), key=lambda r: r["mean_gap"])
        dflt = res[f"K{K}_rho0.25_T0lam/2_T10.02"]
        out[str(K)] = dict(grid=res, best=best["cfg"], best_gap=best["mean_gap"], default_gap=dflt["mean_gap"])
    json.dump(out, open("results/saq_tune.json", "w"), indent=1)

if part == "POOL":
    a_, b_ = map(int, sys.argv[2].split("-")); K = int(sys.argv[3]); warm(K)
    cfg = json.load(open("results/saq_tune.json"))[str(K)]["best"]; jl = f"results/saq_pool_K{K}.jsonl"
    done = {json.loads(l)["seed"] for l in open(jl)} if os.path.exists(jl) else set()
    for seed in range(a_, b_):
        if seed in done: continue
        W = make_instance(200, seed, K); csr = csr_from_W(W); rec = dict(K=K, seed=seed, cfg=cfg, runs={}); k = 0
        saq_run(W, csr, K, cfg, 10, 1)
        for sweeps, Tmax in ((4000, 30), (16000, 100), (64000, 100), (256000, 100)):
            pool = []; need = None
            while need is None or len(pool) < need:
                t0 = time.perf_counter(); f, phi = saq_run(W, csr, K, cfg, sweeps, 34000 + 100000 * seed + k); k += 1
                pool.append(dict(f=f, phi=phi, t=time.perf_counter() - t0))
                if need is None: need = max(4, min(150, 2 * int(Tmax // pool[0]["t"])))
            rec["runs"][str(sweeps)] = pool
        append_jsonl(jl, rec); print("POOL", K, seed, {s_: len(v) for s_, v in rec["runs"].items()}, flush=True)

"""Tuned path-integral SQA (sqa.py) at matched wall-clock time on the instances of Table 7 (200 cells, K = 3 and 4, seeds 0-4).
 TUNE  grid on held-out seeds 100 and 101 at each K: Trotter slices P in {8, 16}, temperature T in {0.05, 0.15, 0.4}, initial transverse
       field G0 in {1, 3, 8} (final 0.01), penalty scale rho in {0.12, 0.25}; 4,000 sweeps, two runs per seed, scored by the repaired
       cost against the better of two native SA runs of 64,000 sweeps (as for the re-tuning of bSB in round 7).  -> results/sqa_tune.json
 POOL  run pools with the tuned setting at 1,000, 4,000 and 16,000 sweeps, sized as for the heated and constraint-preserving bSB pools
       (twice the number of runs that fit in 30, 100 and 100 s), on the machine of Table 7 under a two-process load.
       -> results/sqa_pool_K<K>.jsonl
       TUNE2 extends the grid once, because the best TUNE setting lay on its edge (lowest temperature, field and penalty scale):
       P in {16, 32}, T in {0.02, 0.05}, G0 in {0.3, 1}, rho in {0.06, 0.12}; same seeds and scoring; the best of both grids is used.
Usage: python3 bench_sqa_tuned.py TUNE | TUNE2 | POOL a-b K"""
import sys, os, json, time, itertools, numpy as np
from exp_common import *
from sqa import sqa_solve

part = sys.argv[1]
GRID = list(itertools.product((8, 16), (0.05, 0.15, 0.4), (1.0, 3.0, 8.0), (0.12, 0.25)))


def sqa_run(W, csr, K, cfg, sweeps, rng, seed):
    a, J, lam = build_qubo(W, K, cfg["rho"])
    x, e = sqa_solve(a, J, cfg["P"], cfg["T"], cfg["G0"], 0.01, sweeps, rng, seed)
    f, phi = decode_eval(x, K, csr)[:2]
    return float(f), float(phi)


if part in ("TUNE", "TUNE2"):
    if part == "TUNE2": GRID = list(itertools.product((16, 32), (0.02, 0.05), (0.3, 1.0), (0.06, 0.12)))
    prev = json.load(open("results/sqa_tune.json")) if part == "TUNE2" else {}
    warm(3); warm(4); out = {}
    for K in (3, 4):
        ref = {}
        for s in (100, 101):
            W = make_instance(200, s, K); csr = csr_from_W(W)
            ref[s] = min(sa_native(np.random.default_rng(s + k).integers(0, K, 200), K, *csr, 2.0, 0.02, 64000, 900 + k)[1] for k in range(2))
        res = {}
        for P, T, G0, rho in GRID:
            cfg = dict(P=P, T=T, G0=G0, rho=rho); v = []; ph = []; ts = []
            for s in (100, 101):
                W = make_instance(200, s, K); csr = csr_from_W(W); rng = np.random.default_rng(31000 + s)
                for k in range(2):
                    t0 = time.perf_counter(); f, phi = sqa_run(W, csr, K, cfg, 4000, rng, 500 + k); ts.append(time.perf_counter() - t0)
                    v.append((f / ref[s] - 1) * 100); ph.append(phi)
            key = f"K{K}_P{P}_T{T}_G{G0}_rho{rho}"; res[key] = dict(cfg=cfg, mean_gap=float(np.mean(v)), phi=float(np.mean(ph)), t=float(np.mean(ts)))
            print("TUNE", key, round(float(np.mean(v)), 2), round(float(np.mean(ph)), 2), round(float(np.mean(ts)), 2), flush=True)
        if part == "TUNE2": res = {**prev[str(K)]["grid"], **res}
        best = min(res.values(), key=lambda r: r["mean_gap"])
        dflt = res[f"K{K}_P16_T0.15_G1.0_rho0.25"]
        out[str(K)] = dict(grid=res, best=best["cfg"], best_gap=best["mean_gap"], default_gap=dflt["mean_gap"])
    json.dump(out, open("results/sqa_tune.json", "w"), indent=1)

if part == "POOL":
    a_, b_ = map(int, sys.argv[2].split("-")); K = int(sys.argv[3]); warm(K)
    cfg = json.load(open("results/sqa_tune.json"))[str(K)]["best"]; jl = f"results/sqa_pool_K{K}.jsonl"
    done = {json.loads(l)["seed"] for l in open(jl)} if os.path.exists(jl) else set()
    for seed in range(a_, b_):
        if seed in done: continue
        W = make_instance(200, seed, K); csr = csr_from_W(W); rng = np.random.default_rng(32000 + 10 * seed + K); rec = dict(K=K, seed=seed, cfg=cfg, runs={}); k = 0
        for sweeps, Tmax in ((1000, 30), (4000, 100), (16000, 100)):
            pool = []; need = None
            while need is None or len(pool) < need:
                t0 = time.perf_counter(); f, phi = sqa_run(W, csr, K, cfg, sweeps, rng, 7000 + k); k += 1
                pool.append(dict(f=f, phi=phi, t=time.perf_counter() - t0))
                if need is None: need = max(4, min(150, 2 * int(Tmax // pool[0]["t"])))
            rec["runs"][str(sweeps)] = pool
        append_jsonl(jl, rec); print("POOL", K, seed, {s_: len(v) for s_, v in rec["runs"].items()}, flush=True)

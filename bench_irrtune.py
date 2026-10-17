"""Tune on held-out instances of the irregular layout itself (jitter 0.3, K = 3, 200 cells, seeds 100 and 101, disjoint from
the 30 test instances 0-29), with the scoring used on the regular layout (repaired cost against the better of two frozen native SA runs of
64,000 sweeps, two runs per seed), every solver over its own grid:
 SAQ   bit-flip SA: penalty scale {0.015, 0.03, 0.06, 0.12, 0.25} x T0 {0.15, 0.5, 1, 2, lambda/2} x T1 {0.02, 0.005}, 64,000 sweeps
 SQA   SQA: P {16, 32} x T {0.02, 0.05, 0.15} x G0 {0.3, 1, 3} x penalty scale {0.06, 0.12, 0.25}, 4,000 sweeps
 BSB   bSB: penalty scale {0.015, 0.03, 0.06, 0.12, 0.25}, default dynamics, 12,000 steps (two batches of eight replicas)
 SA    native SA: starting temperature {0.15, 0.3, 0.6, 1.0, 2.0}, T1 = 0.02, 64,000 sweeps
 POOL a-b which   pools on the 30 test instances with the irregular-tuned setting (only for solvers whose setting changed)
        -> results/irrtune_<which>.json, results/irrpool_<which>_K3.jsonl
Usage: python3 bench_irrtune.py SAQ|SQA|BSB|SA | POOL a-b saq|sqa|bsb"""
import sys, os, json, time, itertools, numpy as np
from exp_common import *
from bench_jitter import make_layout
from bench_lowpen import saq_run, sqa_run, bsb_run

part = sys.argv[1]; K = 3; HO = (100, 101)
inst = lambda s: make_layout(200, s, K, 0.3)


def ref_of():
    ref = {}
    for s in HO:
        W = inst(s); csr = csr_from_W(W)
        ref[s] = min(sa_native(np.random.default_rng(s + k).integers(0, K, 200), K, *csr, 2.0, 0.02, 64000, 900 + k)[1] for k in range(2))
    return ref


def tune(name, grid, run1):
    warm(K); ref = ref_of(); res = {}
    for cfg in grid:
        v = []; ph = []
        for s in HO:
            W = inst(s); csr = csr_from_W(W); rng = np.random.default_rng(41000 + s)
            for k in range(2):
                f, phi = run1(W, csr, cfg, rng, 42000 + 10 * s + k); v.append((f / ref[s] - 1) * 100); ph.append(phi)
        key = json.dumps(cfg, sort_keys=True); res[key] = dict(cfg=cfg, mean_gap=float(np.mean(v)), phi=float(np.mean(ph)))
        print(name, key, round(float(np.mean(v)), 2), round(float(np.mean(ph)), 3), flush=True)
    best = min(res.values(), key=lambda r: r["mean_gap"])
    json.dump(dict(grid=res, best=best["cfg"], best_gap=best["mean_gap"]), open(f"results/irrtune_{name.lower()}.json", "w"), indent=1)


if part == "SAQ":
    tune("SAQ", [dict(rho=r, T0=t0, T1=t1) for r, t0, t1 in itertools.product((0.015, 0.03, 0.06, 0.12, 0.25), (0.15, 0.5, 1.0, 2.0, "lam/2"), (0.02, 0.005))],
         lambda W, csr, cfg, rng, sd: saq_run(W, csr, K, cfg, 64000, sd))
if part == "SQA":
    tune("SQA", [dict(P=P, T=T, G0=G, rho=r) for P, T, G, r in itertools.product((16, 32), (0.02, 0.05, 0.15), (0.3, 1.0, 3.0), (0.06, 0.12, 0.25))],
         lambda W, csr, cfg, rng, sd: sqa_run(W, csr, K, cfg, 4000, rng, sd))
if part == "BSB":
    tune("BSB", [dict(rho=r) for r in (0.015, 0.03, 0.06, 0.12, 0.25)], lambda W, csr, cfg, rng, sd: bsb_run(W, csr, K, cfg["rho"], 4, rng))
if part == "SA":
    tune("SA", [dict(T0=t) for t in (0.15, 0.3, 0.6, 1.0, 2.0)],
         lambda W, csr, cfg, rng, sd: (float(sa_native(rng.integers(0, K, 200), K, *csr, cfg["T0"], 0.02, 64000, sd)[1]), 1.0))

if part == "POOL":
    a_, b_ = map(int, sys.argv[2].split("-")); which = sys.argv[3]; warm(K)
    cfg = json.load(open(f"results/irrtune_{which}.json"))["best"]; jl = f"results/irrpool_{which}_K3.jsonl"
    done = {json.loads(l)["seed"] for l in open(jl)} if os.path.exists(jl) else set()
    if len(sys.argv) > 4 and sys.argv[4] == "redo": done = set()   # rerun under two-process load; the analysis keeps the last record
    plan = {"sqa": ((1000, 30), (4000, 100), (16000, 100)), "saq": ((4000, 30), (16000, 100), (64000, 100), (256000, 100)),
            "bsb": ((1, 10), (4, 30), (16, 100))}[which]
    for seed in range(a_, b_):
        if seed in done: continue
        W = inst(seed); csr = csr_from_W(W); rng = np.random.default_rng(43000 + seed); rec = dict(K=K, seed=seed, cfg=cfg, runs={}); k = 0
        for eff, Tmax in plan:
            pool = []; need = None
            while need is None or len(pool) < need:
                t0 = time.perf_counter()
                if which == "sqa": f, phi = sqa_run(W, csr, K, cfg, eff, rng, 44000 + 1000 * seed + k)
                elif which == "saq": f, phi = saq_run(W, csr, K, cfg, eff, 45000 + 100000 * seed + k)
                else: f, phi = bsb_run(W, csr, K, cfg["rho"], eff, rng)
                k += 1; pool.append(dict(f=f, phi=phi, t=time.perf_counter() - t0))
                if need is None: need = max(4, min(150, 2 * int(Tmax // pool[0]["t"])))
            rec["runs"][str(eff)] = pool
        append_jsonl(jl, rec); print("POOL", which, seed, {s_: len(v) for s_, v in rec["runs"].items()}, flush=True)

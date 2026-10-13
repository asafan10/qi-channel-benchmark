"""Top up the K = 4 native SA pools of bench_fresh67.py at 64,000 and 256,000 sweeps, whose first-run time estimate
left them a few runs short of 2 * floor(100 s / mean time per run) on some instances, so that 100 s has an eligible effort.
Extra runs are appended to results/fresh67_topup.jsonl and merged by analyse_fresh67.py.
Usage: python3 topup_fresh67.py SA-tuned|SA-frozen"""
import sys, json, time, numpy as np
from exp_common import *
from bench_fresh67 import one_run
kind = sys.argv[1]; K = 4; n = 200; warm(K)
done = set()
import os
if os.path.exists("results/fresh67_topup.jsonl"):
    for l in open("results/fresh67_topup.jsonl"): r = json.loads(l); done.add((r["kind"], r["seed"]))
base = {}
for f in ("results/fresh67_K4_200a.jsonl", "results/fresh67_K4_200b.jsonl"):
    for l in open(f): r = json.loads(l); base[r["seed"]] = r
for seed in range(5):
    if (kind, seed) in done: continue
    W, csr, a, J, lam, jj = setup(n, seed, K); rng = np.random.default_rng(68000 + seed + (0 if kind == "SA-tuned" else 50)); add = {}
    for e in ("64", "256"):
        pool = list(base[seed]["runs"][kind][e]); new = []
        while len(pool) + len(new) < 2 * int(100 // np.mean([x["t"] for x in pool + new])) + 2:
            new.append(one_run(kind, int(e), n, K, W, csr, a, J, jj, rng, 0))
        add[e] = new
    append_jsonl("results/fresh67_topup.jsonl", dict(K=K, n=n, seed=seed, kind=kind, runs=add))
    print(kind, seed, {e: len(v) for e, v in add.items()}, flush=True)

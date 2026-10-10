"""Fresh, larger run pools for protocol P6 on the irregular layout (jitter 0.3, 200 cells, K = 3, seeds 0-29).
The expected best-of-N estimator is used only where the pool m of stored runs satisfies m >= 2N at the budget, so the pools are sized for the
1-100 s budgets. All runs are new and timed on one machine under the same load (two single-threaded processes on two cores); they are not mixed
with the runs of bench_sigma.py. Settings are those of bench_hybrid.py / bench_hybrid_ext.py / bench_sigma.py, unchanged.
Pools per instance:
  bSB x1 + polish (4k, 16k, 64k sweeps, T 0.3 -> 0.02, all from the same decoded plan): 78 batches
  bSB x4: 48 batches; bSB x16: 40 batches (bSB alone)
  native SA, tuned (T0 = 0.6) and frozen (T0 = 2.0): x1: 50, x4: 40, x16: 34, x64: 26, x256: 20 runs
Usage: python3 bench_sigma_ext.py a-b    -> results/sigma03_ext.jsonl"""
import sys, time
import numpy as np
import exp_common
from bench_jitter import make_layout
exp_common.make_instance = lambda n, seed, K: make_layout(n, seed, K, 0.3)
from exp_common import *

K, n = 3, 200; OUT = "results/sigma03_ext.jsonl"
POL = [4000, 16000, 64000]
BSB_PLAN = [(1, 78), (4, 48), (16, 40)]
SA_PLAN = [(1, 50), (4, 40), (16, 34), (64, 26), (256, 20)]


def run(seed):
    W, csr, a, J, lam, jj = setup(n, seed, K)
    rng = np.random.default_rng(77000 + seed); out = {"bSB": {}, "SA-tuned": {}, "SA-frozen": {}}
    for f, nb in BSB_PLAN:
        recs = []
        for b in range(nb):
            t0 = time.perf_counter(); Xb = sb_batch(a, J, f * CFG["sb_steps"], CFG["sb_R"], rng, False)
            es = [qubo_energy(x, a, *jj) for x in Xb]; sel = int(np.argmin(es))
            ch, _ = repair(Xb[sel].astype(np.int64), K, *csr); fsel = interference(ch, *csr); t_sb = time.perf_counter() - t0
            r = dict(f_sel=float(fsel), t=t_sb, pol={})
            if f == 1:
                for s in POL:
                    t0 = time.perf_counter(); v = sa_native(ch, K, *csr, 0.3, 0.02, s, 1700000 + 1000 * b + s % 997 + 100000 * seed)[1]
                    r["pol"][str(s)] = dict(f=float(v), t=time.perf_counter() - t0)
            recs.append(r)
        out["bSB"][str(f)] = recs
    for key, T0, base in (("SA-tuned", 0.6, 2100000), ("SA-frozen", 2.0, 2900000)):
        for f, nr in SA_PLAN:
            recs = []
            for k in range(nr):
                t0 = time.perf_counter(); v = sa_native(rng.integers(0, K, n), K, *csr, T0, 0.02, f * CFG["native_sweeps"], base + 1000 * f + k + 10000 * seed)[1]
                recs.append(dict(f=float(v), t=time.perf_counter() - t0))
            out[key][str(f)] = recs
    return out


if __name__ == "__main__":
    a_, b_ = map(int, sys.argv[1].split("-")); warm(K); done = done_keys(OUT)
    for seed in range(a_, b_):
        if (n, seed, K) in done: continue
        t0 = time.perf_counter(); res = run(seed)
        append_jsonl(OUT, dict(n=n, seed=seed, K=K, sigma=0.3, runs=res))
        print("S", seed, "%.0fs" % (time.perf_counter() - t0), flush=True)

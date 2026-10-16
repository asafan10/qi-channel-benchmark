"""Full refinement protocol P6 on the irregular layout (jitter 0.3, 200 cells) at K = 4, seeds 0-29: the K = 3 protocol of
bench_sigma_ext.py + bench_sigma_ext2.py (pools of both, combined) repeated at four channels with every setting transferred unchanged
(polish T 0.3 -> 0.02; tuned native SA T0 = 0.15, the value tuned at K = 4 on the main layouts; frozen T0 = 2.0). Two single-threaded
processes on two cores, no other jobs. Decision rule fixed in PREREG_irregular_K4.md before any run.
Pools per instance: bSB x1 + polish (4k, 16k, 64k sweeps) 104 batches; bSB x4 48; bSB x16 56; native SA tuned x1 50, x4 40, x16 34,
x64 26, x256 30; frozen the same with x256 28.
Usage: python3 bench_sigma_k4.py a-b    -> results/sigma03_K4.jsonl"""
import sys, time
import numpy as np
import exp_common
from bench_jitter import make_layout
exp_common.make_instance = lambda n, seed, K: make_layout(n, seed, K, 0.3)
from exp_common import *

K, n = 4, 200; OUT = "results/sigma03_K4.jsonl"
POL = [4000, 16000, 64000]
BSB_PLAN = [(1, 104), (4, 48), (16, 56)]
SA_PLAN = {"SA-tuned": [(1, 50), (4, 40), (16, 34), (64, 26), (256, 30)], "SA-frozen": [(1, 50), (4, 40), (16, 34), (64, 26), (256, 28)]}


def run(seed):
    W, csr, a, J, lam, jj = setup(n, seed, K)
    rng = np.random.default_rng(87000 + seed); out = {"bSB": {}, "SA-tuned": {}, "SA-frozen": {}}
    for f, nb in BSB_PLAN:
        recs = []
        for b in range(nb):
            t0 = time.perf_counter(); Xb = sb_batch(a, J, f * CFG["sb_steps"], CFG["sb_R"], rng, False)
            es = [qubo_energy(x, a, *jj) for x in Xb]; sel = int(np.argmin(es))
            ch, nf = repair(Xb[sel].astype(np.int64), K, *csr); fsel = interference(ch, *csr); t_sb = time.perf_counter() - t0
            r = dict(f_sel=float(fsel), phi=float(decode_eval(Xb[sel], K, csr)[1]), t=t_sb, pol={})
            if f == 1:
                for s in POL:
                    t0 = time.perf_counter(); v = sa_native(ch, K, *csr, 0.3, 0.02, s, 3700000 + 1000 * b + s % 997 + 100000 * seed)[1]
                    r["pol"][str(s)] = dict(f=float(v), t=time.perf_counter() - t0)
            recs.append(r)
        out["bSB"][str(f)] = recs
    for key, T0, base in (("SA-tuned", 0.15, 4100000), ("SA-frozen", 2.0, 4900000)):
        for f, nr in SA_PLAN[key]:
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

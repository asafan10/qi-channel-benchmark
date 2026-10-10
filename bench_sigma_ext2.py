"""Extension of the fresh pools of bench_sigma_ext.py so that 100 s can be tested (analysis plan, amendment 3).
Same instances, settings and machine; two single-threaded processes and no other jobs. Appended runs per instance:
native SA x256 tuned +10 and frozen +8; bSB x16 +16 batches; bSB x1 with 4k/16k/64k polish +26 batches.
Usage: python3 bench_sigma_ext2.py a-b   -> results/sigma03_ext2.jsonl"""
import sys, time
import numpy as np
import exp_common
from bench_jitter import make_layout
exp_common.make_instance = lambda n, seed, K: make_layout(n, seed, K, 0.3)
from exp_common import *

K, n = 3, 200; OUT = "results/sigma03_ext2.jsonl"; POL = [4000, 16000, 64000]


def run(seed):
    W, csr, a, J, lam, jj = setup(n, seed, K)
    rng = np.random.default_rng(99000 + seed); out = {"bSB": {}, "SA-tuned": {}, "SA-frozen": {}}
    for f, nb in ((1, 26), (16, 16)):
        recs = []
        for b in range(nb):
            t0 = time.perf_counter(); Xb = sb_batch(a, J, f * CFG["sb_steps"], CFG["sb_R"], rng, False)
            es = [qubo_energy(x, a, *jj) for x in Xb]; sel = int(np.argmin(es))
            ch, _ = repair(Xb[sel].astype(np.int64), K, *csr); fsel = interference(ch, *csr); t_sb = time.perf_counter() - t0
            r = dict(f_sel=float(fsel), t=t_sb, pol={})
            if f == 1:
                for s in POL:
                    t0 = time.perf_counter(); v = sa_native(ch, K, *csr, 0.3, 0.02, s, 5700000 + 1000 * b + s % 997 + 100000 * seed)[1]
                    r["pol"][str(s)] = dict(f=float(v), t=time.perf_counter() - t0)
            recs.append(r)
        out["bSB"][str(f)] = recs
    for key, T0, base, nr in (("SA-tuned", 0.6, 6100000, 10), ("SA-frozen", 2.0, 6900000, 8)):
        recs = []
        for k in range(nr):
            t0 = time.perf_counter(); v = sa_native(rng.integers(0, K, n), K, *csr, T0, 0.02, 256 * CFG["native_sweeps"], base + k + 10000 * seed)[1]
            recs.append(dict(f=float(v), t=time.perf_counter() - t0))
        out[key]["256"] = recs
    return out


if __name__ == "__main__":
    a_, b_ = map(int, sys.argv[1].split("-")); warm(K); done = done_keys(OUT)
    for seed in range(a_, b_):
        if (n, seed, K) in done: continue
        t0 = time.perf_counter(); res = run(seed)
        append_jsonl(OUT, dict(n=n, seed=seed, K=K, sigma=0.3, runs=res))
        print("S", seed, "%.0fs" % (time.perf_counter() - t0), flush=True)

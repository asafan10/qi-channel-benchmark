"""Full refinement protocol (P6) on the irregular layout of protocol P11, jitter sigma = 0.3.
200 cells, K = 3, seeds 0-29 of bench_jitter.make_layout(n, seed, K, 0.3). Per instance, exactly as for the main layouts:
bench_hybrid.run (bSB x1/x4/x16, descent, polish 1k/4k/16k; native SA x1-x256 frozen), bench_hybrid_ext.run (64k polish),
tuned native SA (T0 = 0.6, the value tuned on the main layouts), the reuse pattern and reuse + native SA polish (reference),
and the dual-corrected semidefinite bound (SCS tolerance 1e-4). Settings tuned on the main layouts are transferred unchanged.
Usage: python3 bench_sigma.py a-b   -> results/sigma03.jsonl"""
import sys, time
import numpy as np
import exp_common
from bench_jitter import make_layout
exp_common.make_instance = lambda n, seed, K: make_layout(n, seed, K, 0.3)      # setup() now builds the sigma = 0.3 layout
from exp_common import *
import bench_hybrid, bench_hybrid_ext
from bench_reuse_polish import best_reuse
from ortools.sat.python import cp_model  # noqa: F401  (import before cvxpy)
from bench_sdp_rig import rig

K, n = 3, 200; OUT = "results/sigma03.jsonl"
if __name__ == "__main__":
    a_, b_ = map(int, sys.argv[1].split("-")); warm(K); done = done_keys(OUT)
    for seed in range(a_, b_):
        if (n, seed, K) in done: continue
        t0 = time.perf_counter(); W, csr, a, J, lam, jj = setup(n, seed, K)
        hyb = bench_hybrid.run(n, seed, K); ext = bench_hybrid_ext.run(n, seed, K)
        rng = np.random.default_rng(43000 + seed); tuned = {}
        for f, nr in bench_hybrid.SA_PLAN:
            recs = []
            for k in range(nr):
                t1 = time.perf_counter(); v = sa_native(rng.integers(0, K, n), K, *csr, 0.6, 0.02, f * CFG["native_sweeps"], 910000 + 1000 * f + k)[1]
                recs.append(dict(f=float(v), t=time.perf_counter() - t1))
            tuned[str(f)] = recs
        p = best_reuse(n, K, csr); f_re = float(interference(p, *csr))
        f_rs = min(float(sa_native(p, K, *csr, T0, 0.02, 64000, 900 + k)[1]) for k, T0 in enumerate((0.3, 0.15)))
        sd = rig(W, K, eps=1e-4)
        append_jsonl(OUT, dict(n=n, seed=seed, K=K, sigma=0.3, hybrid=hyb, ext=ext, tuned=tuned, reuse=f_re, reuse_sa=f_rs, lb=sd["lb_f_rig"], sdp_t=sd["t"]))
        print("S", seed, "%.0fs" % (time.perf_counter() - t0), flush=True)

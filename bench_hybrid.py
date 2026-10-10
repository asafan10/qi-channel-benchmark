"""
bench_hybrid.py -- Experiment A: does native-space refinement after bSB decoding remove the plateau?

Per instance (evaluation seeds 0-4), on one machine so that all timings are comparable:
  bSB     : steps = f*3000, f in {1, 4, 16}, batch of 8 replicas; best-energy replica decoded + repaired (as in the paper)
            + 'desc'  : steepest 1-opt descent in native space on the selected replica
            + 'pol_s' : native SA polish from the selected replica, T0=0.3 -> T1=0.02, s sweeps, s in {1000, 4000, 16000}
              (T0 chosen on held-out seeds 100-101, pilot_polish.py; frozen)
  SA-native: sweeps = f*1000, f in {1, 4, 16, 64, 256}, random start, T0=2.0 -> 0.02 (frozen baseline)
Each record stores cost and wall-clock of every component so that matched-wall-clock analysis can add them.
"""
import sys, time
import numpy as np
from exp_common import *

BSB_PLAN = [(1, 6), (4, 4), (16, 2)]
POL = [1000, 4000, 16000]; POL_T0, POL_T1 = 0.3, 0.02
SA_PLAN = [(1, 40), (4, 20), (16, 10), (64, 5), (256, 3)]


def run(n, seed, K):
    W, csr, a, J, lam, jj = setup(n, seed, K)
    rng = np.random.default_rng(42000 + seed); out = {"bSB": {}, "SA": {}}
    for f, nb in BSB_PLAN:
        recs = []
        for b in range(nb):
            t0 = time.perf_counter(); Xb = sb_batch(a, J, f * CFG["sb_steps"], CFG["sb_R"], rng, False)
            es = [qubo_energy(x, a, *jj) for x in Xb]; sel = int(np.argmin(es))
            ch, _ = repair(Xb[sel].astype(np.int64), K, *csr); fsel = interference(ch, *csr); t_sb = time.perf_counter() - t0
            t0 = time.perf_counter(); _, fd, mv = descent(ch, K, *csr); t_d = time.perf_counter() - t0
            r = dict(f_sel=float(fsel), t=t_sb, f_desc=float(fd), t_desc=t_d, moves=int(mv), pol={})
            for s in POL:
                t0 = time.perf_counter(); v = sa_native(ch, K, *csr, POL_T0, POL_T1, s, 800000 + 1000 * b + s % 997)[1]
                r["pol"][str(s)] = dict(f=float(v), t=time.perf_counter() - t0)
            recs.append(r)
        out["bSB"][str(f)] = recs
    for f, nr in SA_PLAN:
        recs = []
        for k in range(nr):
            t0 = time.perf_counter(); v = sa_native(rng.integers(0, K, n), K, *csr, 2.0, 0.02, f * CFG["native_sweeps"], 900000 + 1000 * f + k)[1]
            recs.append(dict(f=float(v), t=time.perf_counter() - t0))
        out["SA"][str(f)] = recs
    return out


if __name__ == "__main__":
    SEEDS = range(*map(int, sys.argv[4].split("-"))) if len(sys.argv) > 4 else range(5)
    K = int(sys.argv[1]); sizes = [int(x) for x in sys.argv[2].split(",")]; outf = sys.argv[3]
    warm(K); done = done_keys(outf)
    for n in sizes:
        for seed in SEEDS:
            if (n, seed, K) in done: continue
            t0 = time.perf_counter(); res = run(n, seed, K)
            best = min([r["f_sel"] for v in res["bSB"].values() for r in v] + [p["f"] for v in res["bSB"].values() for r in v for p in r["pol"].values()]
                       + [r["f"] for v in res["SA"].values() for r in v])
            append_jsonl(outf, dict(n=n, seed=seed, K=K, best=best, runs=res))
            print(K, n, seed, "%.0fs" % (time.perf_counter() - t0), "best %.2f" % best, flush=True)

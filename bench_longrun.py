"""
bench_longrun.py -- symmetric long-run reference (Priority 1).

bSB      : steps = f * 3000, f in {4, 8, 16, 32, 64}  (12k ... 192k steps), batch of 8 replicas
SA-native: sweeps = f * 1000, f in {64, 256}            (64k, 256k sweeps)
Per run we store decoded costs of all 8 bSB replicas (and the best-energy replica's cost) plus wall-clock.
"""
import json, sys, os, time
import numpy as np
from bench_ran_qubo import *
from bench_tts import warmup, previous_bks

PLAN = {"bSB": [(4, 5), (8, 3), (16, 3), (32, 2), (64, 1)],
        "SA-native": [(64, 8), (256, 3)]}


def run_instance(n, seed, K, cfg):
    W = make_instance(n, seed, K); csr = csr_from_W(W)
    a, J, lam = build_qubo(W, K, cfg["lam_scale"])
    jp, jx, jd = J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64)
    rng = np.random.default_rng(31000 + seed); out = {"bSB": {}, "SA-native": {}}
    for f, nr in PLAN["bSB"]:
        recs = []
        for k in range(nr):
            t0 = time.perf_counter()
            Xb = sb_batch(a, J, int(f * cfg["sb_steps"]), cfg["sb_R"], rng, False)
            es = [qubo_energy(x, a, jp, jx, jd) for x in Xb]
            fr = [float(decode_eval(x, K, csr)[0]) for x in Xb]
            recs.append(dict(f_sel=fr[int(np.argmin(es))], f_all=fr, t=time.perf_counter() - t0))
        out["bSB"][str(f)] = recs
    for f, nr in PLAN["SA-native"]:
        recs = []
        for k in range(nr):
            t0 = time.perf_counter()
            v = sa_native(rng.integers(0, K, n), K, *csr, 2.0, 0.02, int(f * cfg["native_sweeps"]), 700000 + 1000 * f + k)[1]
            recs.append(dict(f_sel=float(v), t=time.perf_counter() - t0))
        out["SA-native"][str(f)] = recs
    return out


if __name__ == "__main__":
    K = int(sys.argv[1]); sizes = [int(x) for x in sys.argv[2].split(",")]; seeds = int(sys.argv[3]); outf = sys.argv[4]
    warmup(K, CFG); prev = previous_bks(K); done = set()
    if os.path.exists(outf):
        for l in open(outf):
            d = json.loads(l); done.add((d["n"], d["seed"]))
    with open(outf, "a") as fh:
        for n in sizes:
            for s in range(seeds):
                if (n, s) in done: continue
                t0 = time.perf_counter(); res = run_instance(n, s, K, CFG)
                best = min(min(r["f_sel"] for r in v) for sv in res.values() for v in sv.values())
                best = min(best, min(min(r["f_all"]) for v in res["bSB"].values() for r in v))
                fh.write(json.dumps(dict(n=n, seed=s, K=K, prev_bks=prev[(n, s)], new_best=best, runs=res)) + "\n"); fh.flush()
                print(n, s, "%.0fs" % (time.perf_counter() - t0), "prev_bks %.2f longrun_best %.2f" % (prev[(n, s)], best), flush=True)

"""Fresh, one-machine run pools for the long-run and matched wall-clock comparison (Tables 6 and 7, Figures 1-2).
Native SA with the frozen (T0 = 2) and the tuned starting temperature (T0 = 0.6 at K = 3, 0.15 at K = 4; T1 = 0.02), untuned bSB
(dt = 0.5, linear pump, 8 replicas, best-energy replica decoded) and, at K = 3, the frozen tuned bSB configuration of Section 5.3.
Each effort e gets a pool of runs large enough for the exact expected best of N without replacement up to its largest budget
Tmax(e): at least 2 * floor(Tmax / t_e) runs (and at least 4), t_e measured on the first run. Solvers are interleaved run by run
within an instance, so any load from the second process affects every solver alike.
Usage: python3 bench_fresh67.py K n seed_from-seed_to out.jsonl"""
import sys, time, math, numpy as np
from exp_common import *
from bench_round5 import sb2

T0_TUNED = {3: 0.6, 4: 0.15}
SA_PLAN = [(4, 10), (16, 30), (64, 100), (256, 100)]          # (multiple of 1,000 sweeps, Tmax in s)
SB_PLAN = [(1, 10), (4, 30), (16, 100), (64, 100)]           # (multiple of 3,000 steps, Tmax in s)
SBT_PLAN = [(4, 30), (16, 100)]
TUNED_CFG = dict(dt=0.5, gamma=2.0, c0m=2.0)


def one_run(kind, f, n, K, W, csr, a, J, jj, rng, k):
    t0 = time.perf_counter()
    if kind in ("SA-frozen", "SA-tuned"):
        T0 = 2.0 if kind == "SA-frozen" else T0_TUNED[K]
        v = sa_native(rng.integers(0, K, n), K, *csr, T0, 0.02, f * CFG["native_sweeps"], int(rng.integers(1, 2**31 - 1)))[1]
        return dict(f=float(v), t=time.perf_counter() - t0)
    if kind == "bSB":
        Xb = sb_batch(a, J, f * CFG["sb_steps"], CFG["sb_R"], rng, False)
    else:
        Xb = sb2(a, J, f * CFG["sb_steps"], CFG["sb_R"], rng, TUNED_CFG["dt"], TUNED_CFG["gamma"], TUNED_CFG["c0m"])
    es = [qubo_energy(x, a, *jj) for x in Xb]; x = Xb[int(np.argmin(es))]
    v, phi = decode_eval(x, K, csr)[:2]
    return dict(f=float(v), t=time.perf_counter() - t0)


if __name__ == "__main__":
    K = int(sys.argv[1]); n = int(sys.argv[2]); a_, b_ = map(int, sys.argv[3].split("-")); outf = sys.argv[4]
    warm(K); done = done_keys(outf)
    plans = {"SA-frozen": SA_PLAN, "SA-tuned": SA_PLAN, "bSB": SB_PLAN}
    if K == 3: plans["bSB-tuned"] = SBT_PLAN
    for seed in range(a_, b_):
        if (n, seed, K) in done: continue
        W, csr, a, J, lam, jj = setup(n, seed, K); rng = np.random.default_rng(67000 + 100 * K + seed)
        pools = {(kd, f): [] for kd, pl in plans.items() for f, _ in pl}
        need = {(kd, f): None for kd, pl in plans.items() for f, _ in pl}
        tmax = {(kd, f): T for kd, pl in plans.items() for f, T in pl}
        t_start = time.time()
        while True:
            open_ = [key for key in pools if need[key] is None or len(pools[key]) < need[key]]
            if not open_: break
            for key in open_:                       # one run of every unfinished effort per round (interleaved)
                r = one_run(key[0], key[1], n, K, W, csr, a, J, jj, rng, len(pools[key])); pools[key].append(r)
                if need[key] is None: need[key] = max(4, 2 * int(tmax[key] // r["t"]))
        out = {kd: {str(f): pools[(kd, f)] for f, _ in pl} for kd, pl in plans.items()}
        append_jsonl(outf, dict(n=n, seed=seed, K=K, wall=time.time() - t_start, runs=out))
        print(K, n, seed, "done in", round(time.time() - t_start), "s", {f"{k[0]}x{k[1]}": len(v) for k, v in pools.items()}, flush=True)

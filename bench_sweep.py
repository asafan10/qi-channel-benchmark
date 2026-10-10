"""
bench_sweep.py -- run-length sweep for independent-run TTS.

For each solver the per-run effort is varied over a grid (multiples of the frozen baseline):
  Tabu      : iterations   = f * 40 n          f in {0.25, 1, 4, 16}
  SA-native : sweeps       = f * 1000          f in {0.25, 1, 4, 16}
  SA-QUBO   : sweeps       = f * 3000          f in {0.25, 1, 4, 16}
  PT-QUBO   : sweeps/replica = f * 400 (12 replicas)  f in {0.25, 1, 4, 16}
  bSB / dSB : steps        = f * 3000 (batch of 8)    f in {1/6, 1/3, 1, 2}  (500, 1000, 3000, 6000 steps)
The number of independent runs per cell adapts to a wall-clock budget (min 8, max 300).
"""
import json, sys, os, time
import numpy as np
from bench_ran_qubo import *
from bench_tts import warmup, previous_bks

SOLVERS = ["Tabu", "SA-native", "SA-QUBO", "PT-QUBO", "bSB", "dSB"]
FACT = {"Tabu": [0.25, 1, 4, 16], "SA-native": [0.25, 1, 4, 16], "SA-QUBO": [0.25, 1, 4, 16],
        "PT-QUBO": [0.25, 1, 4, 16], "bSB": [1 / 6, 1 / 3, 1, 2], "dSB": [1 / 6, 1 / 3, 1, 2]}
TBUDGET = {"Tabu": 5.0, "SA-native": 5.0, "bSB": 5.0, "SA-QUBO": 3.0, "PT-QUBO": 3.0, "dSB": 3.0}
RMIN, RMAX = 8, 300


def one_run(s, f, n, K, csr, a, J, jp, jx, jd, lam, temps, rng, sd, cfg):
    N = n * K; T0, T1 = cfg["T0_scale"] * lam, cfg["T1"]
    if s == "Tabu":
        return tabu(rng.integers(0, K, n), K, *csr, max(1, int(round(f * cfg["tabu_iters_per_cell"] * n))), 8, sd)[1]
    if s == "SA-native":
        return sa_native(rng.integers(0, K, n), K, *csr, 2.0, 0.02, max(1, int(round(f * cfg["native_sweeps"]))), sd)[1]
    if s == "SA-QUBO":
        return decode_eval(sa_qubo(rng.integers(0, 2, N), a, jp, jx, jd, T0, T1, max(2, int(round(f * cfg["sweeps"]))), sd)[0], K, csr)[0]
    if s == "PT-QUBO":
        return decode_eval(pt_qubo(rng.integers(0, 2, (cfg["pt_R"], N)), a, jp, jx, jd, temps, max(2, int(round(f * cfg["pt_sweeps"]))), sd)[0], K, csr)[0]
    Xb = sb_batch(a, J, int(round(f * cfg["sb_steps"])), cfg["sb_R"], rng, s == "dSB")
    es = [qubo_energy(x, a, jp, jx, jd) for x in Xb]
    return decode_eval(Xb[int(np.argmin(es))], K, csr)[0]


def run_instance(n, seed, K, cfg):
    W = make_instance(n, seed, K); csr = csr_from_W(W)
    a, J, lam = build_qubo(W, K, cfg["lam_scale"])
    jp, jx, jd = J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64)
    temps = np.geomspace(cfg["pt_tmin"], cfg["T0_scale"] * lam, cfg["pt_R"])
    rng = np.random.default_rng(9000 + seed); out = {}
    for si, s in enumerate(SOLVERS):
        out[s] = {}
        for fi, f in enumerate(FACT[s]):
            base = 10 ** 6 * seed + 10 ** 4 * si + 100 * fi
            fs, ts = [], []
            k = 0
            while True:
                t0 = time.perf_counter()
                v = one_run(s, f, n, K, csr, a, J, jp, jx, jd, lam, temps, rng, base + k, cfg)
                ts.append(time.perf_counter() - t0); fs.append(float(v)); k += 1
                if k >= RMAX: break
                if k >= RMIN and sum(ts) >= TBUDGET[s]: break
            out[s][f"{f:.4f}"] = dict(f=fs, t=ts)
    return out


if __name__ == "__main__":
    K = int(sys.argv[1]); sizes = [int(x) for x in sys.argv[2].split(",")]; seeds = int(sys.argv[3]); outf = sys.argv[4]
    warmup(K, CFG); prev = previous_bks(K)
    done = set()
    if os.path.exists(outf):
        for l in open(outf):
            d = json.loads(l); done.add((d["n"], d["seed"]))
    with open(outf, "a") as fh:
        for n in sizes:
            for s in range(seeds):
                if (n, s) in done: continue
                t0 = time.perf_counter()
                res = run_instance(n, s, K, CFG)
                best_new = min(min(c["f"]) for sv in res.values() for c in sv.values())
                fh.write(json.dumps(dict(n=n, seed=s, K=K, prev_bks=prev[(n, s)], new_best=best_new, runs=res)) + "\n"); fh.flush()
                print(n, s, "%.0fs" % (time.perf_counter() - t0), "prev_bks %.2f sweep_best %.2f" % (prev[(n, s)], best_new), flush=True)

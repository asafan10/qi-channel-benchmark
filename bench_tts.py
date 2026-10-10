"""
bench_tts.py -- independent-run protocol for time-to-solution (TTS) with bootstrap CIs.

For each instance (same generator/seeds as bench_ran_qubo.py) and solver we perform R
INDEPENDENT runs. A "run" is one call of the solver's configured routine:
  Tabu      : one tabu search, 40n iterations
  SA-native : one anneal, 1000 sweeps
  SA-QUBO   : one anneal on the QUBO, 3000 sweeps
  PT-QUBO   : one parallel-tempering run, 12 replicas x 400 sweeps
  bSB / dSB : one batch of 8 bifurcation replicas, 3000 steps (best-energy replica decoded)
Per run we store the decoded cost f and wall-clock time t. Hyper-parameters = frozen CFG.
"""
import json, sys, os, time
import numpy as np
from bench_ran_qubo import *

SOLVERS = ["Tabu", "SA-native", "SA-QUBO", "PT-QUBO", "bSB", "dSB"]
RUNS = {"Tabu": 300, "SA-native": 200, "bSB": 30, "SA-QUBO": 40, "PT-QUBO": 20, "dSB": 15}


def previous_bks(K):
    files = {3: ["results_K3.jsonl", "results_K3_fix.jsonl"], 4: ["results_K4.jsonl"]}[K]
    m = {}
    for f in files:
        for l in open(f):
            r = json.loads(l); key = (r["n"], r["seed"])
            m.setdefault(key, {}).update({s: v["E"] for s, v in r["res"].items()})
    return {k: min(v.values()) for k, v in m.items()}


def warmup(K, cfg):
    W = make_instance(12, 0, K); csr = csr_from_W(W); a, J, lam = build_qubo(W, K, cfg["lam_scale"])
    jp, jx, jd = J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64)
    rng = np.random.default_rng(0); N = 12 * K
    tabu(rng.integers(0, K, 12), K, *csr, 50, 8, 0); sa_native(rng.integers(0, K, 12), K, *csr, 2.0, 0.02, 5, 0)
    x = sa_qubo(rng.integers(0, 2, N), a, jp, jx, jd, 1.0, 0.1, 5, 0)[0]; decode_eval(x, K, csr)
    pt_qubo(rng.integers(0, 2, (3, N)), a, jp, jx, jd, np.geomspace(0.1, 1.0, 3), 3, 0)
    Xb = sb_batch(a, J, 20, 2, rng, False); qubo_energy(Xb[0], a, jp, jx, jd)


def run_instance(n, seed, K, cfg):
    W = make_instance(n, seed, K); csr = csr_from_W(W)
    a, J, lam = build_qubo(W, K, cfg["lam_scale"])
    jp, jx, jd = J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64)
    N = n * K; rng = np.random.default_rng(5000 + seed)
    T0, T1 = cfg["T0_scale"] * lam, cfg["T1"]
    temps = np.geomspace(cfg["pt_tmin"], T0, cfg["pt_R"])
    out = {}
    for si, s in enumerate(SOLVERS):
        fs, ts = [], []
        for k in range(RUNS[s]):
            sd = 100000 * seed + 1000 * si + k
            t0 = time.perf_counter()
            if s == "Tabu":
                f = tabu(rng.integers(0, K, n), K, *csr, cfg["tabu_iters_per_cell"] * n, 8, sd)[1]
            elif s == "SA-native":
                f = sa_native(rng.integers(0, K, n), K, *csr, 2.0, 0.02, cfg["native_sweeps"], sd)[1]
            elif s == "SA-QUBO":
                f = decode_eval(sa_qubo(rng.integers(0, 2, N), a, jp, jx, jd, T0, T1, cfg["sweeps"], sd)[0], K, csr)[0]
            elif s == "PT-QUBO":
                f = decode_eval(pt_qubo(rng.integers(0, 2, (cfg["pt_R"], N)), a, jp, jx, jd, temps, cfg["pt_sweeps"], sd)[0], K, csr)[0]
            else:
                Xb = sb_batch(a, J, cfg["sb_steps"], cfg["sb_R"], rng, s == "dSB")
                es = [qubo_energy(x, a, jp, jx, jd) for x in Xb]
                f = decode_eval(Xb[int(np.argmin(es))], K, csr)[0]
            ts.append(time.perf_counter() - t0); fs.append(float(f))
        out[s] = dict(f=fs, t=ts)
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
                res = run_instance(n, s, K, CFG)
                best_new = min(min(v["f"]) for v in res.values())
                fh.write(json.dumps(dict(n=n, seed=s, K=K, prev_bks=prev[(n, s)], new_best=best_new, runs=res)) + "\n"); fh.flush()
                print(n, s, "prev_bks %.2f new_best %.2f" % (prev[(n, s)], best_new),
                      {k: round(float(np.mean(v["t"])), 3) for k, v in res.items()}, flush=True)

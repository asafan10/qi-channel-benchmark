"""Cheap follow-ups to bench_round6.py.
 K2B  dual-corrected semidefinite lower bounds (bench_sdp_rig.rig, SCS tolerance 1e-4) for the K = 2 instances of part K2, so that the
      two-channel gaps can be stated against a bound and not only against the best plan found.
 MG   matched-gap control for the phase-domain analysis: native SA plans at 1,000, 2,000 and 4,000 sweeps (four runs, seeds 0-4, regular
      layout, K = 3), whose gap is close to that of bSB, are stored for analyse_round6.py.
 TUNE the penalty scales of the domain-wall encoding and of the per-cell penalty re-chosen on three held-out instances (seeds
      100-102) instead of one.
 TIME per-run time of native SA (64,000 sweeps, 20 runs per process) and bSB (12,000 steps, 10 batches per process) with one process and
      with two concurrent processes on the two cores of machine B, on three instances (K = 3 at 200 and 500 cells, K = 4 at 200; run alone).
Usage: python3 bench_round6b.py K2B|MG|TUNE|TIME -> results/round6b_<part>.json"""
import sys, json, time, subprocess, numpy as np
import exp_common
from exp_common import *
from bench_jitter import make_layout

part = sys.argv[1]; n = 200; out = {}


def layout(seed, K, sig):
    return make_instance(n, seed, K) if sig == 0.15 else make_layout(n, seed, K, sig)


if part == "K2B":
    from bench_sdp_rig import rig
    rows = []
    for sig in (0.15, 0.3):
        for seed in range(5):
            W = layout(seed, 2, sig); r = rig(W, 2, eps=1e-4)
            rows.append(dict(sigma=sig, seed=seed, lb=r["lb_f_rig"], t=r["t"], status=r["status"]))
            print("K2B", sig, seed, round(r["lb_f_rig"], 3), round(r["t"], 1), flush=True)
    out["rows"] = rows
if part == "MG":
    warm(3); rows = []
    for seed in range(5):
        W = layout(seed, 3, 0.15); csr = csr_from_W(W); rng = np.random.default_rng(6500 + seed); rec = dict(seed=seed, plans={})
        for f in (1, 2, 4):
            rec["plans"][f"sa{f}"] = [np.asarray(sa_native(rng.integers(0, 3, n), 3, *csr, 2.0, 0.02, f * 1000, int(rng.integers(1, 2**31 - 1)))[0]).tolist()
                                      for _ in range(4)]
        rows.append(rec); print("MG", seed, flush=True)
    out["rows"] = rows
if part == "TUNE":
    import importlib.util
    spec = importlib.util.spec_from_file_location("r6", "bench_round6.py"); src = open("bench_round6.py").read()
    g = {"__name__": "r6"}; exec(src[:src.index("part = sys.argv[1]")], g)
    res = {"dw": {}, "pc": {}}
    for K in (3, 4):
        sc = {}
        for rho in (0.25, 0.5, 1.0, 2.0):
            v = []
            for s in (100, 101, 102):
                W = g["layout"](s, K, 0.15); csr = csr_from_W(W); a, J, lam = g["dw_qubo"](W, K, rho)
                X, _ = g["best_of_batches"](a, J, 3000, np.random.default_rng(10 + s))
                v.append(np.mean([interference(g["dw_decode"](x, K)[0], *csr) for x in X]) / g["sa_runs"](W, K, 64, np.random.default_rng(s), nr=2)[0][0][1])
            sc[rho] = float(np.mean(v))
        res["dw"][str(K)] = sc; print("TUNE dw", K, sc, flush=True)
    sc = {}
    for rho in (0.25, 0.5, 1.0):
        v = []
        for s in (100, 101, 102):
            W = g["layout"](s, 3, 0.15); csr = csr_from_W(W)
            lam_i = rho * (W.sum(1) + 1.0); N = n * 3; a = np.repeat(-lam_i, 3); quad = {}
            ii, jj = np.nonzero(np.triu(W, 1))
            for i, j in zip(ii, jj):
                for c in range(3): quad[(i * 3 + c, j * 3 + c)] = W[i, j]
            for i in range(n):
                for c in range(3):
                    for cp_ in range(c + 1, 3): quad[(i * 3 + c, i * 3 + cp_)] = 2 * lam_i[i]
            a, J = g["qubo_from_terms"](N, {p: a[p] for p in range(N)}, quad)
            X, _ = g["best_of_batches"](a, J, 3000, np.random.default_rng(10 + s))
            v.append(np.mean([interference(repair(x.astype(np.int64), 3, *csr)[0], *csr) for x in X]) / g["sa_runs"](W, 3, 64, np.random.default_rng(s), nr=2)[0][0][1])
        sc[rho] = float(np.mean(v))
    res["pc"]["3"] = sc; print("TUNE pc", sc, flush=True); out = res
if part == "TIME":
    # worker: python3 bench_round6b.py TIME K n seed rngseed -> per-run times of 20 native SA runs (64,000 sweeps) and 10 bSB batches (12,000 steps)
    if len(sys.argv) > 3:
        K_, n_, s_ = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); warm(K_); W, csr, a, J, lam, jj = setup(n_, s_, K_)
        rng = np.random.default_rng(int(sys.argv[5])); ts = {"sa": [], "bsb": []}
        for _ in range(20):
            t0 = time.perf_counter(); sa_native(rng.integers(0, K_, n_), K_, *csr, 2.0, 0.02, 64000, int(rng.integers(1, 2**31 - 1))); ts["sa"].append(time.perf_counter() - t0)
        for _ in range(10):
            t0 = time.perf_counter(); sb_batch(a, J, 12000, 8, rng, False); ts["bsb"].append(time.perf_counter() - t0)
        print(json.dumps(ts)); sys.exit(0)
    def run(k, inst):
        ps = [subprocess.Popen([sys.executable, "bench_round6b.py", "TIME", *map(str, inst), str(100 + i)], stdout=subprocess.PIPE, text=True) for i in range(k)]
        return [json.loads(p.communicate()[0].strip().splitlines()[-1]) for p in ps]
    def summ(L, key):
        x = np.array([v for d in L for v in d[key]]); return dict(n_runs=int(len(x)), mean=float(x.mean()), sd=float(x.std(ddof=1)), cv=float(x.std(ddof=1) / x.mean()))
    out = {}
    for inst in ((3, 200, 0), (3, 500, 0), (4, 200, 0)):
        one = run(1, inst); two = run(2, inst)
        out["%d_%d_%d" % inst] = {"one": {k: summ(one, k) for k in ("sa", "bsb")}, "two": {k: summ(two, k) for k in ("sa", "bsb")}}
        print(inst, json.dumps(out["%d_%d_%d" % inst]), flush=True)
json.dump(out, open(f"results/round6b_{part}.json", "w"), indent=1)

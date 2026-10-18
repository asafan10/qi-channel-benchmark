"""Irregular layout (jitter 0.3, K = 3, 200 cells; test seeds 0-29, held-out seeds 100-101):
 PTTUNE   parallel tempering on the QUBO (PT-QUBO) tuned on the held-out instances with the scoring of bench_irrtune.py:
          replicas {8, 16} x penalty scale {0.03, 0.06, 0.12} x lowest temperature {0.02, 0.05, 0.1} x highest {0.5, 1.0, lambda/2},
          2,000 sweeps, one-hot random starts, best visited state decoded; only settings with raw feasibility >= 0.95 qualify
 PTTUNE2  extends that grid where its best setting lay on the edge (replicas {16, 32} x penalty scale {0.06, 0.09} x lowest
          temperature {0.1, 0.2, 0.3} x highest {0.5, 1.0})
 PTTUNE3  extends it once more (replicas {32, 64} x penalty scale {0.09, 0.12, 0.18} x lowest temperature {0.3, 0.4, 0.6} x highest
          {0.5, 1, 2}, highest above lowest)
 PTPOOL a-b   pools of the tuned PT-QUBO on the test instances (500, 2,000 and 8,000 sweeps; pools sized for 10, 30 and 30 s)
 POLSTART a-b the polish of protocol P6 (native SA, T 0.3 -> 0.02, 16,000 and 64,000 sweeps; four runs each) started from: a random plan,
          the row-order repair of the empty vector, greedy by weighted degree, weighted DSATUR, the best of 1,000 descended roundings of the
          semidefinite relaxation (SCS, tolerance 1e-4), and the decoded state of one bSB batch (a new batch for each run)
 LONGPOL a-b  bSB x1 followed by a polish of 256,000 sweeps (the longest effort of native SA), 24 runs per instance (pools for 30 and 100 s)
-> results/r12_*.json(l)
Usage: python3 bench_round12.py PTTUNE | PTPOOL a-b | POLSTART a-b | LONGPOL a-b"""
import os, sys, json, time, itertools
from ortools.sat.python import cp_model  # noqa: F401  (OR-Tools before CVXPY and Numba)
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from exp_common import *
from bench_jitter import make_layout
from bench_ran_qubo import pt_qubo, greedy
from decoder_floor import order_greedy, dsatur

K = 3; HO = (100, 101); inst = lambda s: make_layout(200, s, K, 0.3)


def onehot(n, rng, R):
    X = np.zeros((R, n * K), dtype=np.int64)
    for r in range(R): X[r, np.arange(n) * K + rng.integers(0, K, n)] = 1
    return X


def pt_run(W, csr, cfg, sweeps, rng, seed):
    a, J, lam = build_qubo(W, K, cfg["rho"]); jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64))
    tmax = 0.5 * lam if cfg["Tmax"] == "lam/2" else float(cfg["Tmax"])
    temps = np.geomspace(cfg["Tmin"], max(tmax, 1.01 * cfg["Tmin"]), cfg["R"])
    x, e = pt_qubo(onehot(W.shape[0], rng, cfg["R"]), a, *jj, temps, sweeps, seed)[:2]
    f, phi = decode_eval(x, K, csr)[:2]
    return float(f), float(phi)


def ref_of():
    ref = {}
    for s in HO:
        W = inst(s); csr = csr_from_W(W)
        ref[s] = min(sa_native(np.random.default_rng(s + k).integers(0, K, 200), K, *csr, 2.0, 0.02, 64000, 900 + k)[1] for k in range(2))
    return ref


part = sys.argv[1]
if part == "PTTUNE":
    warm(K); ref = ref_of(); res = {}
    for R, rho, tmin, tmax in itertools.product((8, 16), (0.03, 0.06, 0.12), (0.02, 0.05, 0.1), (0.5, 1.0, "lam/2")):
        cfg = dict(R=R, rho=rho, Tmin=tmin, Tmax=tmax); v = []; ph = []
        for s in HO:
            W = inst(s); csr = csr_from_W(W); rng = np.random.default_rng(51000 + s)
            for k in range(2):
                f, phi = pt_run(W, csr, cfg, 2000, rng, 52000 + 10 * s + k); v.append((f / ref[s] - 1) * 100); ph.append(phi)
        key = json.dumps(cfg, sort_keys=True); res[key] = dict(cfg=cfg, mean_gap=float(np.mean(v)), phi=float(np.mean(ph)))
        print("PT", key, round(float(np.mean(v)), 2), round(float(np.mean(ph)), 3), flush=True)
    best = min([r for r in res.values() if r["phi"] >= 0.95] or res.values(), key=lambda r: r["mean_gap"])
    json.dump(dict(grid=res, best=best["cfg"], best_gap=best["mean_gap"]), open("results/r12_pttune.json", "w"), indent=1)

if part == "PTTUNE2":   # the first grid's best lay on its edge (16 replicas, lowest temperature 0.1): extend it there, same scoring
    warm(K); ref = ref_of(); d = json.load(open("results/r12_pttune.json")); res = d["grid"]
    for R, rho, tmin, tmax in itertools.product((16, 32), (0.06, 0.09), (0.1, 0.2, 0.3), (0.5, 1.0)):
        cfg = dict(R=R, rho=rho, Tmin=tmin, Tmax=tmax); key = json.dumps(cfg, sort_keys=True)
        if key in res: continue
        v = []; ph = []
        for s in HO:
            W = inst(s); csr = csr_from_W(W); rng = np.random.default_rng(51000 + s)
            for k in range(2):
                f, phi = pt_run(W, csr, cfg, 2000, rng, 52000 + 10 * s + k); v.append((f / ref[s] - 1) * 100); ph.append(phi)
        res[key] = dict(cfg=cfg, mean_gap=float(np.mean(v)), phi=float(np.mean(ph)))
        print("PT2", key, round(float(np.mean(v)), 2), round(float(np.mean(ph)), 3), flush=True)
    best = min([r for r in res.values() if r["phi"] >= 0.95] or res.values(), key=lambda r: r["mean_gap"])
    json.dump(dict(grid=res, best=best["cfg"], best_gap=best["mean_gap"], extended=True), open("results/r12_pttune.json", "w"), indent=1)

if part == "PTTUNE3":   # the extended grid's best still lay on its edge (32 replicas, lowest temperature 0.3); extend once more
    warm(K); ref = ref_of(); d = json.load(open("results/r12_pttune.json")); res = d["grid"]
    for R, rho, tmin, tmax in itertools.product((32, 64), (0.09, 0.12, 0.18), (0.3, 0.4, 0.6), (0.5, 1.0, 2.0)):
        if tmax <= tmin: continue
        cfg = dict(R=R, rho=rho, Tmin=tmin, Tmax=tmax); key = json.dumps(cfg, sort_keys=True)
        if key in res: continue
        v = []; ph = []
        for s in HO:
            W = inst(s); csr = csr_from_W(W); rng = np.random.default_rng(51000 + s)
            for k in range(2):
                f, phi = pt_run(W, csr, cfg, 2000, rng, 52000 + 10 * s + k); v.append((f / ref[s] - 1) * 100); ph.append(phi)
        res[key] = dict(cfg=cfg, mean_gap=float(np.mean(v)), phi=float(np.mean(ph)))
        print("PT3", key, round(float(np.mean(v)), 2), round(float(np.mean(ph)), 3), flush=True)
    best = min([r for r in res.values() if r["phi"] >= 0.95] or res.values(), key=lambda r: r["mean_gap"])
    json.dump(dict(grid=res, best=best["cfg"], best_gap=best["mean_gap"], extended=True, extended2=True, best_before=d["best"]),
              open("results/r12_pttune.json", "w"), indent=1)

if part == "PTPOOL":
    a_, b_ = map(int, sys.argv[2].split("-")); warm(K); jl = "results/r12_ptpool.jsonl"
    while not os.path.exists("results/r12_pttune.json") or not json.load(open("results/r12_pttune.json")).get("extended"): time.sleep(5)
    cfg = json.load(open("results/r12_pttune.json"))["best"]
    done = {json.loads(l)["seed"] for l in open(jl)} if os.path.exists(jl) else set()
    for seed in range(a_, b_):
        if seed in done: continue
        W = inst(seed); csr = csr_from_W(W); rng = np.random.default_rng(53000 + seed); rec = dict(K=K, seed=seed, cfg=cfg, runs={}); k = 0
        for eff, Tmax in ((500, 10), (2000, 30), (8000, 30)):
            pool = []; need = None
            while need is None or len(pool) < need:
                t0 = time.perf_counter(); f, phi = pt_run(W, csr, cfg, eff, rng, 54000 + 1000 * seed + k); k += 1
                pool.append(dict(f=f, phi=phi, t=time.perf_counter() - t0))
                if need is None: need = max(4, min(150, 2 * int(Tmax // pool[0]["t"])))
            rec["runs"][str(eff)] = pool
        append_jsonl(jl, rec); print("PTPOOL", seed, {s_: len(v) for s_, v in rec["runs"].items()}, flush=True)

if part == "POLSTART":
    from sdp_round import solve_sdp   # imports OR-Tools before CVXPY (library load order)
    a_, b_ = map(int, sys.argv[2].split("-")); warm(K); jl = "results/r12_polstart.jsonl"
    done = {json.loads(l)["seed"] for l in open(jl)} if os.path.exists(jl) else set()
    for seed in range(a_, b_):
        if seed in done: continue
        W = inst(seed); csr = csr_from_W(W); n = W.shape[0]; a, J, lam = build_qubo(W, K, CFG["lam_scale"])
        jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64)); rng = np.random.default_rng(55000 + seed)
        t0 = time.perf_counter(); X, _ = solve_sdp(W, K, eps=1e-4)
        w_, V = np.linalg.eigh((X + X.T) / 2); V = V * np.sqrt(np.clip(w_, 0, None)); best = None
        for r in range(1000):
            g = rng.standard_normal((V.shape[1], K)); ch = np.argmax(V @ g, axis=1).astype(np.int64)
            ch, e, _ = descent(ch, K, *csr)
            if best is None or e < best[1]: best = (ch.copy(), e)
        t_sdp = time.perf_counter() - t0
        starts = {"row-order": order_greedy(np.arange(n), K, *csr), "greedy": greedy(csr, n, K), "DSATUR": dsatur(K, *csr), "SDP rounding": best[0]}
        rec = dict(K=K, seed=seed, t_sdp=t_sdp, start_cost={k_: float(interference(np.asarray(v, dtype=np.int64), *csr)) for k_, v in starts.items()}, runs={})
        for name in ("random", "row-order", "greedy", "DSATUR", "SDP rounding", "bSB"):
            rec["runs"][name] = {}
            for sw in (16000, 64000):
                rr = []
                for k in range(4):
                    t0 = time.perf_counter(); f_start = None
                    if name == "random": ch = rng.integers(0, K, n).astype(np.int64)
                    elif name == "bSB":
                        Xb = sb_batch(a, J, CFG["sb_steps"], CFG["sb_R"], rng, False)
                        x = Xb[int(np.argmin([qubo_energy(x, a, *jj) for x in Xb]))]
                        ch, _ = repair(x.astype(np.int64), K, *csr); f_start = float(interference(ch, *csr))
                    else: ch = np.asarray(starts[name], dtype=np.int64).copy()
                    v = sa_native(ch, K, *csr, 0.3, 0.02, sw, 56000 + 100000 * seed + 1000 * k + sw % 997)[1]
                    rr.append(dict(f=float(v), t=time.perf_counter() - t0, f_start=f_start))
                rec["runs"][name][str(sw)] = rr
        append_jsonl(jl, rec); print("POLSTART", seed, {k_: round(v, 1) for k_, v in rec["start_cost"].items()}, round(t_sdp), flush=True)

if part == "LONGPOL":
    a_, b_ = map(int, sys.argv[2].split("-")); warm(K); jl = "results/r12_longpol.jsonl"
    done = {json.loads(l)["seed"] for l in open(jl)} if os.path.exists(jl) else set()
    if len(sys.argv) > 3 and sys.argv[3] == "redo": done = set()
    for seed in range(a_, b_):
        if seed in done: continue
        W = inst(seed); csr = csr_from_W(W); a, J, lam = build_qubo(W, K, CFG["lam_scale"])
        jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64)); rng = np.random.default_rng(57000 + seed); rr = []
        for k in range(24):
            t0 = time.perf_counter(); Xb = sb_batch(a, J, CFG["sb_steps"], CFG["sb_R"], rng, False)
            x = Xb[int(np.argmin([qubo_energy(x, a, *jj) for x in Xb]))]; ch, _ = repair(x.astype(np.int64), K, *csr)
            v = sa_native(ch, K, *csr, 0.3, 0.02, 256000, 58000 + 100000 * seed + k)[1]
            rr.append(dict(f=float(v), t=time.perf_counter() - t0))
        append_jsonl(jl, dict(K=K, seed=seed, runs=rr)); print("LONGPOL", seed, round(float(np.median([r["f"] for r in rr])), 2), round(float(np.mean([r["t"] for r in rr])), 1), flush=True)

"""Equal tuning of the penalty for bSB and bit-flip SA, and the single-bit annealers in the 30-instance protocol.
 SAQ2   extends the bit-flip SA grid of bench_saq_tuned.py below its boundary optimum: penalty scale 0.03 and 0.015 with the same
        temperatures (held-out seeds 100-101, 64,000 sweeps, same scoring); merged into results/saq_tune.json (the best of both grids is used).
 BSBT   bSB with the penalty scale tuned on the same held-out instances and scoring: rho in {0.015, 0.03, 0.06, 0.12, 0.25}, default
        dynamics (dt 0.5, linear pump, 8 replicas, best-energy replica decoded and repaired), 12,000 steps, two batches per seed.
        -> results/bsb_rho_tune.json
 BSBP a-b K rho   bSB pools at penalty scale rho on the instances of Table 7 (200 cells, seeds 0-4), efforts and pool sizes of
        bench_fresh67.py (1x/4x/16x/64x 3,000 steps for budgets up to 10/30/100/100 s); raw feasibility recorded.
        -> results/bsb_rho<rho>_K<K>.jsonl
 IRR a-b which   tuned SQA or tuned bit-flip SA (settings of results/sqa_tune.json and results/saq_tune.json at K = 3, transferred unchanged)
        on the 30 irregular instances of Table 12 (jitter 0.3, K = 3, 200 cells), pools as in bench_sqa_tuned.py / bench_saq_tuned.py.
        -> results/irr_<which>_K3.jsonl
Every job runs in a two-process load (run_round7.sh)."""
import sys, os, json, time, itertools, numpy as np
from exp_common import *
from sqa import sqa_solve

part = sys.argv[1]


def saq_run(W, csr, K, cfg, sweeps, seed):
    a, J, lam = build_qubo(W, K, cfg["rho"]); jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64))
    T0 = 0.5 * lam if cfg["T0"] == "lam/2" else float(cfg["T0"])
    x0 = np.zeros(len(a), dtype=np.int64); ch = np.random.default_rng(seed).integers(0, K, W.shape[0]); x0[np.arange(W.shape[0]) * K + ch] = 1
    bx, _ = sa_qubo(x0, a, *jj, T0, cfg["T1"], sweeps, seed)
    f, phi = decode_eval(bx, K, csr)[:2]
    return float(f), float(phi)


def sqa_run(W, csr, K, cfg, sweeps, rng, seed):
    a, J, lam = build_qubo(W, K, cfg["rho"])
    x, e = sqa_solve(a, J, cfg["P"], cfg["T"], cfg["G0"], 0.01, sweeps, rng, seed)
    f, phi = decode_eval(x, K, csr)[:2]
    return float(f), float(phi)


def bsb_run(W, csr, K, rho, f, rng):
    a, J, lam = build_qubo(W, K, rho); jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64))
    Xb = sb_batch(a, J, f * CFG["sb_steps"], CFG["sb_R"], rng, False)
    es = [qubo_energy(x, a, *jj) for x in Xb]; x = Xb[int(np.argmin(es))]
    v, phi = decode_eval(x, K, csr)[:2]
    return float(v), float(phi)


def heldout_ref(K):
    ref = {}
    for s in (100, 101):
        W = make_instance(200, s, K); csr = csr_from_W(W)
        ref[s] = min(sa_native(np.random.default_rng(s + k).integers(0, K, 200), K, *csr, 2.0, 0.02, 64000, 900 + k)[1] for k in range(2))
    return ref


if part == "SAQ2":
    warm(3); warm(4); tu = json.load(open("results/saq_tune.json"))
    for K in (3, 4):
        ref = heldout_ref(K); res = {}
        for rho, T0, T1 in itertools.product((0.03, 0.015), (0.15, 0.5, 1.0, 2.0, "lam/2"), (0.02, 0.005)):
            cfg = dict(rho=rho, T0=T0, T1=T1); v = []; ph = []; ts = []
            for s in (100, 101):
                W = make_instance(200, s, K); csr = csr_from_W(W)
                for k in range(2):
                    t0 = time.perf_counter(); f, phi = saq_run(W, csr, K, cfg, 64000, 33000 + 10 * s + k); ts.append(time.perf_counter() - t0)
                    v.append((f / ref[s] - 1) * 100); ph.append(phi)
            key = f"K{K}_rho{rho}_T0{T0}_T1{T1}"; res[key] = dict(cfg=cfg, mean_gap=float(np.mean(v)), phi=float(np.mean(ph)), t=float(np.mean(ts)))
            print("SAQ2", key, round(float(np.mean(v)), 2), round(float(np.mean(ph)), 3), flush=True)
        g = {**tu[str(K)]["grid"], **res}; best = min(g.values(), key=lambda r: r["mean_gap"])
        tu[str(K)].update(grid=g, best=best["cfg"], best_gap=best["mean_gap"], extended=True)
    json.dump(tu, open("results/saq_tune.json", "w"), indent=1)

if part == "BSBT":
    warm(3); warm(4); out = {}
    for K in (3, 4):
        ref = heldout_ref(K); res = {}
        for rho in (0.015, 0.03, 0.06, 0.12, 0.25):
            v = []; ph = []
            for s in (100, 101):
                W = make_instance(200, s, K); csr = csr_from_W(W); rng = np.random.default_rng(35000 + s)
                for k in range(2):
                    f, phi = bsb_run(W, csr, K, rho, 4, rng); v.append((f / ref[s] - 1) * 100); ph.append(phi)
            res[str(rho)] = dict(mean_gap=float(np.mean(v)), phi=float(np.mean(ph)))
            print("BSBT", K, rho, round(float(np.mean(v)), 2), round(float(np.mean(ph)), 3), flush=True)
        out[str(K)] = dict(grid=res, best=float(min(res, key=lambda r: res[r]["mean_gap"])))
    json.dump(out, open("results/bsb_rho_tune.json", "w"), indent=1)

if part == "BSBP":
    a_, b_ = map(int, sys.argv[2].split("-")); K = int(sys.argv[3]); rho = float(sys.argv[4]); warm(K)
    jl = f"results/bsb_rho{rho:g}_K{K}.jsonl"; done = {json.loads(l)["seed"] for l in open(jl)} if os.path.exists(jl) else set()
    for seed in range(a_, b_):
        if seed in done: continue
        W = make_instance(200, seed, K); csr = csr_from_W(W); rng = np.random.default_rng(36000 + 100 * K + seed); rec = dict(K=K, seed=seed, rho=rho, runs={})
        for f, Tmax in ((1, 10), (4, 30), (16, 100), (64, 100)):
            pool = []; need = None
            while need is None or len(pool) < need:
                t0 = time.perf_counter(); v, phi = bsb_run(W, csr, K, rho, f, rng); pool.append(dict(f=v, phi=phi, t=time.perf_counter() - t0))
                if need is None: need = max(4, 2 * int(Tmax // pool[0]["t"]))
            rec["runs"][str(f)] = pool
        append_jsonl(jl, rec); print("BSBP", K, rho, seed, {k: len(v) for k, v in rec["runs"].items()}, flush=True)

if part == "IRR":
    from bench_jitter import make_layout
    a_, b_ = map(int, sys.argv[2].split("-")); which = sys.argv[3]; K = 3; warm(K)
    jl = f"results/irr_{which}_K3.jsonl"; done = {json.loads(l)["seed"] for l in open(jl)} if os.path.exists(jl) else set()
    cfg = json.load(open("results/sqa_tune.json" if which == "sqa" else "results/saq_tune.json"))["3"]["best"]
    plan = ((1000, 30), (4000, 100), (16000, 100)) if which == "sqa" else ((4000, 30), (16000, 100), (64000, 100), (256000, 100))
    for seed in range(a_, b_):
        if seed in done: continue
        W = make_layout(200, seed, K, 0.3); csr = csr_from_W(W); rng = np.random.default_rng(37000 + seed); rec = dict(K=K, seed=seed, cfg=cfg, runs={}); k = 0
        for sweeps, Tmax in plan:
            pool = []; need = None
            while need is None or len(pool) < need:
                t0 = time.perf_counter()
                f, phi = (sqa_run(W, csr, K, cfg, sweeps, rng, 38000 + 1000 * seed + k) if which == "sqa" else saq_run(W, csr, K, cfg, sweeps, 39000 + 100000 * seed + k))
                k += 1; pool.append(dict(f=f, phi=phi, t=time.perf_counter() - t0))
                if need is None: need = max(4, min(150, 2 * int(Tmax // pool[0]["t"])))
            rec["runs"][str(sweeps)] = pool
        append_jsonl(jl, rec); print("IRR", which, seed, {s_: len(v) for s_, v in rec["runs"].items()}, flush=True)

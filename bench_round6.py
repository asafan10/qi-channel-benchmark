"""Three controls for the encoding and the plateau (n = 200; seeds 0-4; regular layout sigma_j = 0.15 and
irregular layout sigma_j = 0.3 of protocol P11).
 K2  K = 2 on the same layouts: one binary variable per cell, no one-hot constraint, zero Ising field (plain weighted Max-Cut).
     bSB (3,000 x {1, 4, 16} steps, 4 batches of 8 replicas, best-energy replica) against native SA (1,000 x {4, 16, 64, 256}
     sweeps, 4 runs, T0 = 2 -> 0.02). Gaps against the best plan found by any run of this part.
 DW  bSB on a domain-wall encoding (Chancellor 2019): K - 1 binary variables per cell, y_0 = 1 and y_K = 0 fixed, channel c
     selected by d_c = y_{c-1} - y_c, one penalty lambda * y_c (1 - y_{c-1}) per interior wall position (K = 3: one per cell),
     lambda = rho_dw (Delta_w + 1) with rho_dw chosen on held-out seed 100 from {0.25, 0.5, 1, 2}. Decoding takes the number of
     leading ones as the channel, so every sample gives a plan without repair. Compared with one-hot bSB (rho = 0.25) and native
     SA at the same steps or sweeps, on both layouts at K = 3 and on the regular layout at K = 4. Gaps against the reference of
     the main experiments (regular) or of Table 12 (irregular).
 PC  one-hot bSB with a per-cell penalty lambda_i = rho (delta_i + 1), rho chosen on held-out seed 100 from {0.25, 0.5, 1}.
 PH  plans of one-hot bSB, DW-bSB and native SA from part DW (K = 3) are stored for the phase-domain analysis (analyse_round6.py).
Usage: python3 bench_round6.py K2|DW|PC -> results/round6_<part>.json"""
import sys, json, math, time, numpy as np, scipy.sparse as sp
import exp_common
from exp_common import *
from bench_jitter import make_layout

SEEDS = range(5); n = 200


def layout(seed, K, sig):
    return make_instance(n, seed, K) if sig == 0.15 else make_layout(n, seed, K, sig)


def best_of_batches(a, J, steps, rng, nb=4):
    jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64)); out = []; ts = []
    for b in range(nb):
        t0 = time.perf_counter(); Xb = sb_batch(a, J, steps, CFG["sb_R"], rng, False)
        es = [qubo_energy(x, a, *jj) for x in Xb]; out.append(Xb[int(np.argmin(es))]); ts.append(time.perf_counter() - t0)
    return out, float(np.mean(ts))


def sa_runs(W, K, f, rng, nr=4):
    csr = csr_from_W(W); res = []; ts = []
    for k in range(nr):
        t0 = time.perf_counter(); ch, v = sa_native(rng.integers(0, K, n), K, *csr, 2.0, 0.02, f * 1000, int(rng.integers(1, 2**31 - 1)))
        res.append((np.asarray(ch).tolist(), float(v))); ts.append(time.perf_counter() - t0)
    return res, float(np.mean(ts))


def qubo_from_terms(N, lin, quad):
    a = np.zeros(N)
    for v, c in lin.items(): a[v] += c
    rows, cols, vals = [], [], []
    for (u, v), c in quad.items(): rows += [u, v]; cols += [v, u]; vals += [c, c]
    return a, sp.csr_matrix((vals, (rows, cols)), shape=(N, N))


def dw_qubo(W, K, rho):
    """Domain-wall QUBO: variable index i*(K-1) + (c-1) for y_{i,c}, c = 1..K-1."""
    m = K - 1; N = n * m; lam = rho * (W.sum(1).max() + 1.0)
    def form(i, c):   # d_{i,c} = y_{c-1} - y_c as (const, {var: coef}), c = 1..K
        const = 0.0; lin = {}
        if c - 1 == 0: const += 1.0
        else: lin[i * m + c - 2] = lin.get(i * m + c - 2, 0) + 1.0
        if c == K: pass
        else: lin[i * m + c - 1] = lin.get(i * m + c - 1, 0) - 1.0
        return const, lin
    lin = {}; quad = {}; E0 = 0.0
    ii, jj = np.nonzero(np.triu(W, 1))
    for i, j in zip(ii, jj):
        w = W[i, j]
        for c in range(1, K + 1):
            (ci, li), (cj, lj) = form(i, c), form(j, c)
            E0 += w * ci * cj
            for v, x in li.items(): lin[v] = lin.get(v, 0) + w * x * cj
            for v, x in lj.items(): lin[v] = lin.get(v, 0) + w * x * ci
            for u, x in li.items():
                for v, y in lj.items(): key = (min(u, v), max(u, v)); quad[key] = quad.get(key, 0) + w * x * y
    for i in range(n):
        for c in range(2, K):     # penalty lam * y_c (1 - y_{c-1})
            u, v = i * m + c - 2, i * m + c - 1
            lin[v] = lin.get(v, 0) + lam; quad[(u, v)] = quad.get((u, v), 0) - lam
    a, J = qubo_from_terms(N, lin, quad); return a, J, lam


def dw_decode(y, K):
    Y = y.reshape(n, K - 1); ch = np.zeros(n, np.int64); valid = np.ones(n, bool)
    for i in range(n):
        k = 0
        while k < K - 1 and Y[i, k] == 1: k += 1
        ch[i] = k; valid[i] = Y[i, k:].sum() == 0
    return ch, float(valid.mean())


def k2_qubo(W):
    lin = {}; quad = {}
    ii, jj = np.nonzero(np.triu(W, 1))
    for i, j in zip(ii, jj):
        w = W[i, j]; lin[i] = lin.get(i, 0) - w; lin[j] = lin.get(j, 0) - w; quad[(i, j)] = 2 * w
    return qubo_from_terms(n, lin, quad)


def irregular_ref():
    src = open("analyse_sigma_ext.py").read(); g = {}
    exec(src[:src.index("LB = {")], g); return g["REF"]


part = sys.argv[1]; out = {}
warm(3); warm(4)
if part == "K2":
    rows = []
    for sig in (0.15, 0.3):
        for seed in SEEDS:
            W = layout(seed, 2, sig); csr = csr_from_W(W); a, J = k2_qubo(W); rng = np.random.default_rng(6200 + seed + int(100 * sig))
            rec = dict(sigma=sig, seed=seed, bsb={}, sa={})
            for f in (1, 4, 16):
                X, t = best_of_batches(a, J, f * 3000, rng)
                rec["bsb"][str(f)] = dict(f=[float(interference(x.astype(np.int64), *csr)) for x in X], t=t)
            for f in (4, 16, 64, 256):
                res, t = sa_runs(W, 2, f, rng); rec["sa"][str(f)] = dict(f=[v for _, v in res], t=t)
            rec["best"] = min(min(v["f"]) for d in (rec["bsb"], rec["sa"]) for v in d.values())
            rows.append(rec); print("K2", sig, seed, round(rec["best"], 3), {k: round(min(v["f"]) / rec["best"] * 100 - 100, 2) for k, v in rec["bsb"].items()},
                                   {k: round(min(v["f"]) / rec["best"] * 100 - 100, 2) for k, v in rec["sa"].items()}, flush=True)
    out["rows"] = rows
if part == "DW":
    from refs import build_ref
    REFR = build_ref(); REFI = irregular_ref()
    tune = {}
    for K in (3, 4):          # held-out choice of rho_dw on seed 100, regular layout, 3,000 steps
        W = layout(100, K, 0.15); csr = csr_from_W(W); sc = {}
        for rho in (0.25, 0.5, 1.0, 2.0):
            a, J, lam = dw_qubo(W, K, rho); X, _ = best_of_batches(a, J, 3000, np.random.default_rng(1))
            sc[rho] = float(np.mean([interference(dw_decode(x, K)[0], *csr) for x in X]))
        tune[K] = min(sc, key=sc.get); print("tune", K, sc, "->", tune[K], flush=True)
    rows = []
    for K, sig in ((3, 0.15), (3, 0.3), (4, 0.15)):
        for seed in SEEDS:
            W = layout(seed, K, sig); csr = csr_from_W(W); rng = np.random.default_rng(6300 + 10 * K + seed + int(100 * sig))
            ref = REFR[(K, n, seed)] if sig == 0.15 else REFI[seed]
            rec = dict(K=K, sigma=sig, seed=seed, ref=float(ref), dw={}, oh={}, sa={}, plans={})
            a, J, lam = dw_qubo(W, K, tune[K])
            for f in (1, 4, 16):
                X, t = best_of_batches(a, J, f * 3000, rng); dec = [dw_decode(x, K) for x in X]
                rec["dw"][str(f)] = dict(f=[float(interference(c, *csr)) for c, _ in dec], valid=[v for _, v in dec], t=t)
                if f in (1, 16): rec["plans"][f"dw{f}"] = [c.tolist() for c, _ in dec]
            a1, J1, _ = build_qubo(W, K, CFG["lam_scale"])
            for f in (1, 4, 16):
                X, t = best_of_batches(a1, J1, f * 3000, rng); dec = [repair(x.astype(np.int64), K, *csr) for x in X]
                rec["oh"][str(f)] = dict(f=[float(interference(c, *csr)) for c, _ in dec], t=t)
                if f in (1, 16): rec["plans"][f"oh{f}"] = [np.asarray(c).tolist() for c, _ in dec]
            for f in (4, 16, 64, 256):
                res, t = sa_runs(W, K, f, rng); rec["sa"][str(f)] = dict(f=[v for _, v in res], t=t)
                if f in (64, 256): rec["plans"][f"sa{f}"] = [c for c, _ in res]
            g = lambda v: round((np.median(v) / ref - 1) * 100, 2)
            rows.append(rec); print("DW", K, sig, seed, "dw", {k: g(v["f"]) for k, v in rec["dw"].items()}, "valid", round(np.mean(rec["dw"]["1"]["valid"]), 3),
                                   "oh", {k: g(v["f"]) for k, v in rec["oh"].items()}, "sa", {k: g(v["f"]) for k, v in rec["sa"].items()}, flush=True)
    out["tune"] = {str(k): v for k, v in tune.items()}; out["rows"] = rows
if part == "PC":           # per-cell penalty lambda_i = rho (delta_i + 1) instead of the global rho (Delta_w + 1)
    from refs import build_ref
    REFR = build_ref(); REFI = irregular_ref()
    def build_qubo_pc(W, K, rho):
        lam_i = rho * (W.sum(1) + 1.0); N = n * K; a = np.repeat(-lam_i, K); quad = {}
        ii, jj = np.nonzero(np.triu(W, 1))
        for i, j in zip(ii, jj):
            for c in range(K): quad[(i * K + c, j * K + c)] = W[i, j]
        for i in range(n):
            for c in range(K):
                for cp in range(c + 1, K): quad[(i * K + c, i * K + cp)] = 2 * lam_i[i]
        return qubo_from_terms(N, {p: a[p] for p in range(N)}, quad)
    W = layout(100, 3, 0.15); csr = csr_from_W(W); sc = {}
    for rho in (0.25, 0.5, 1.0):
        a, J = build_qubo_pc(W, 3, rho); X, _ = best_of_batches(a, J, 3000, np.random.default_rng(1))
        sc[rho] = float(np.mean([interference(repair(x.astype(np.int64), 3, *csr)[0], *csr) for x in X]))
    rho_pc = min(sc, key=sc.get); print("tune pc", sc, "->", rho_pc, flush=True)
    rows = []
    for sig in (0.15, 0.3):
        for seed in SEEDS:
            W = layout(seed, 3, sig); csr = csr_from_W(W); rng = np.random.default_rng(6400 + seed + int(100 * sig))
            ref = REFR[(3, n, seed)] if sig == 0.15 else REFI[seed]; a, J = build_qubo_pc(W, 3, rho_pc); rec = dict(sigma=sig, seed=seed, ref=float(ref), pc={})
            for f in (1, 4, 16):
                X, t = best_of_batches(a, J, f * 3000, rng); dec = [repair(x.astype(np.int64), 3, *csr) for x in X]
                rec["pc"][str(f)] = dict(f=[float(interference(c, *csr)) for c, _ in dec], phi=[nf / n for _, nf in dec], t=t)
            rows.append(rec); print("PC", sig, seed, {k: round((np.median(v["f"]) / ref - 1) * 100, 2) for k, v in rec["pc"].items()}, flush=True)
    out["rho_pc"] = rho_pc; out["tune_pc"] = {str(k): v for k, v in sc.items()}; out["rows"] = rows
json.dump(out, open(f"results/round6_{part}.json", "w"))

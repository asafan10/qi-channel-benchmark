"""Cheap probes of the bSB dynamics at 200 cells (irregular layout, jitter 0.3, K = 3, rho = 0.25).
 A  bSB started from the semidefinite rounding (seeds 0-9): the best descended rounding (sdp_round.py pipeline) is loaded into the
    oscillators (z = 0.9 s, y = 0) and bSB runs 3,000 steps with the pump started at a0 in {0, 0.5, 0.9}; we record whether the
    decoded plan keeps, improves or loses the starting plan.
 B  three hypotheses of Section 6.2 (seeds 0-4, 3,000 steps, best of each batch, 4 batches): heated bSB (momentum heating as in the
    open-source package, coefficient 0.06), a separate field scale c_h / c0 in {0.5, 2} (the force becomes -(c0 J s + c_h h)), and the
    number of replicas R in {8, 32, 128}.
 C  label-order dependence of the repair (seeds 0-4, 100 and 200 cells, regular layout, rho = 0.06 where raw feasibility is below 1):
    SA on the QUBO and bSB outputs repaired in index order and in 5 random cell orders.
 D  replica diversity on G22 (Table B6): mean pairwise Hamming distance between the 32 replicas and number of distinct cuts at
    1,000-30,000 steps.
Usage: python3 bench_round3.py A|B|C|D  -> results/round3_<part>.json"""
import sys, json, math, time, numpy as np, scipy.sparse as sp
import exp_common
from exp_common import *
from bench_ran_qubo import sb_batch, qubo_energy, sa_qubo, build_qubo, make_instance, csr_from_W
from bench_jitter import make_layout


def sb_custom(a, J, steps, R, rng, c_h=1.0, heat=0.0, X0=None, a0=0.0):
    N = len(a); Jm = J * 0.25; h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
    c0 = 0.5 * math.sqrt(N - 1) / math.sqrt((Jm.multiply(Jm)).sum())
    X = rng.uniform(-0.1, 0.1, (N, R)) if X0 is None else np.repeat(X0[:, None], R, 1).astype(float)
    Y = rng.uniform(-0.1, 0.1, (N, R)) if X0 is None else np.zeros((N, R))
    for t in range(steps):
        at = a0 + (1 - a0) * t / steps; Yp = Y.copy()
        Y += 0.5 * (-(1.0 - at) * X - (c0 * (Jm @ X) + c0 * c_h * h[:, None]))
        if heat: Y += 0.5 * heat * Yp
        X += 0.5 * Y
        m = np.abs(X) > 1; X[m] = np.sign(X[m]); Y[m] = 0.0
    return (X > 0).astype(np.int64).T


def decode_best(X, a, jj, K, csr):
    es = [qubo_energy(x, a, *jj) for x in X]; x = X[int(np.argmin(es))]
    ch, _ = repair(x.astype(np.int64), K, *csr); return float(interference(ch, *csr)), float(np.mean(x.reshape(-1, K).sum(1) == 1))


part = sys.argv[1]; out = {}
if part in ("A", "B"):
    exp_common.make_instance = lambda n, seed, K: make_layout(n, seed, K, 0.3)
    REF = json.load(open("results/sigma03_ext.json"))
    OLD = {r["seed"]: r for r in map(json.loads, open("results/sigma03.jsonl"))}
if part == "A":
    from sdp_round import solve_sdp
    rows = []
    for seed in range(10):
        W, csr, a, J, lam, jj = setup(200, seed, 3); K = 3
        Xv, t_sdp = solve_sdp(W, K, 1e-4); w_, U = np.linalg.eigh((Xv + Xv.T) / 2); V = U * np.sqrt(np.clip(w_, 0, None))
        rng = np.random.default_rng(1000 + seed); best = None
        for _ in range(1000):
            ch = np.argmax(V @ rng.normal(size=(200, K)), axis=1).astype(np.int64); chd, fd, _m = descent(ch, K, *csr)
            if best is None or fd < best[1]: best = (chd, fd)
        x0 = np.zeros(200 * K); x0[np.arange(200) * K + best[0]] = 1; z0 = 0.9 * (2 * x0 - 1)
        rec = dict(seed=seed, f_start=float(best[1]), t_sdp=t_sdp, runs={})
        for a0 in (0.0, 0.5, 0.9):
            fs = [decode_best(sb_custom(a, J, 3000, 8, np.random.default_rng(seed * 10 + b), X0=z0, a0=a0), a, jj, K, csr)[0] for b in range(4)]
            rec["runs"][str(a0)] = fs
        rows.append(rec); print("A", seed, rec["f_start"], {k: [round(v, 2) for v in vs] for k, vs in rec["runs"].items()}, flush=True)
    out["rows"] = rows
if part == "B":
    rows = []
    for seed in range(5):
        W, csr, a, J, lam, jj = setup(200, seed, 3); K = 3; rec = dict(seed=seed, res={})
        for name, kw, R in (("default", {}, 8), ("heated", dict(heat=0.06), 8), ("c_h=0.5", dict(c_h=0.5), 8), ("c_h=2", dict(c_h=2.0), 8),
                            ("R=32", {}, 32), ("R=128", {}, 128)):
            vals = []; t0 = time.perf_counter()
            for b in range(4): vals.append(decode_best(sb_custom(a, J, 3000, R, np.random.default_rng(seed * 100 + b), **kw), a, jj, K, csr))
            rec["res"][name] = dict(f=[v[0] for v in vals], phi=[v[1] for v in vals], t=(time.perf_counter() - t0) / 4)
        rows.append(rec); print("B", seed, {k: round(min(v["f"]), 2) for k, v in rec["res"].items()}, flush=True)
    out["rows"] = rows
if part == "C":
    rows = []
    for n in (100, 200):
        for seed in range(5):
            K = 3; W = make_instance(n, seed, K); a, J, lam = build_qubo(W, K, 0.06)
            jp, jx, jd = J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64); csr = csr_from_W(W)
            rng = np.random.default_rng(seed); outs = {}
            xs = {"SA-QUBO": sa_qubo(rng.integers(0, 2, n * K).astype(np.int64), a, jp, jx, jd, 0.5 * lam, 0.15, 3000, seed)[0],
                  "bSB": sb_batch(a, J, 3000, 8, rng, False)[0]}
            for name, x in xs.items():
                x = np.asarray(x, dtype=np.int64); ch, nf = repair(x, K, *csr); f_idx = interference(ch, *csr); f_perm = []
                for p in range(5):
                    perm = np.random.default_rng(100 + p).permutation(n); inv = np.argsort(perm)
                    Wp = W[np.ix_(perm, perm)]; csrp = csr_from_W(Wp); xp = x.reshape(n, K)[perm].ravel()
                    chp, _ = repair(xp, K, *csrp); f_perm.append(float(interference(chp[inv], *csr)))
                outs[name] = dict(phi=float(nf / n), f_index=float(f_idx), f_perm=f_perm)
            rows.append(dict(n=n, seed=seed, res=outs)); print("C", n, seed, {k: (round(v["phi"], 2), round(v["f_index"], 2), round(float(np.mean(v["f_perm"])), 2)) for k, v in outs.items()}, flush=True)
    out["rows"] = rows
if part == "D":
    from gset_check import read_gset, cut
    Wm = read_gset("gset/G22"); n = Wm.shape[0]; deg = np.asarray(Wm.sum(1)).ravel(); a = -deg; J = (2.0 * Wm).tocsr(); rng = np.random.default_rng(5)
    for steps in (1000, 3000, 10000, 30000):
        X = np.vstack([sb_batch(a, J, steps, 8, rng, False) for _ in range(4)]); S = 2.0 * X - 1.0
        cuts = [float(cut(Wm, s)) for s in S]; D = [(S[i] != S[j]).mean() for i in range(len(S)) for j in range(i + 1, len(S))]
        # a cut and its global flip are the same partition: count the Hamming distance up to the flip
        D = [min(d_, 1 - d_) for d_ in D]
        out[str(steps)] = dict(best=max(cuts), median=float(np.median(cuts)), distinct_cuts=len(set(cuts)), mean_hamming=float(np.mean(D)))
        print("D", steps, out[str(steps)], flush=True)
json.dump(out, open(f"results/round3_{part}.json", "w"), indent=1)

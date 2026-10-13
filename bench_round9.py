"""Cheap controls of the bSB pipeline and an exploratory constraint-preserving bSB. 200 cells, main layouts.
 REP   decoding choices (K = 3, 4; seeds 0-4; 3,000, 12,000 and 48,000 steps; four batches of eight replicas): the plan of the
       best-energy replica's final state (the default), the best repaired plan over all eight final states, and the best plan over the
       states visited along each replica's trajectory (sign state recorded every 50 steps; best by energy and best by repaired cost),
       which matches native SA, whose sa_native returns the best state it visited.
 OH2   one-hot bSB at K = 2 (seeds 0-4, rho = 0.25 as in the main experiments): the field h_p = lambda (K - 2) / 2 + delta_i / 4 loses its
       penalty part but keeps the one-hot constraint; compared with the one-variable Max-Cut encoding of round6_K2 on the same layouts.
 POLT  time per sweep of native SA, frozen schedule (T 2 -> 0.02) from a random plan against the polish schedule (T 0.3 -> 0.02) from a
       decoded bSB plan, 64,000 sweeps, seeds 0-4, K = 3; acceptance rate measured by a Python twin over the first and last 1,000 sweeps.
 CPT   constraint-preserving bSB, tuning: per cell the K positions live on the probability simplex; the pump acts around its centre 1/K,
       the force is the interference gradient sum_j w_ij x_jc (no penalty, no field), projected on the simplex's tangent space; after the
       position update each cell is projected back on the simplex (Euclidean projection) and the momentum of every component at the
       boundary x = 0 is set to zero (the analogue of bSB's inelastic walls). Plan = argmax per cell. Grid: time step {0.25, 0.5, 1} x
       coupling multiplier xi {0.5, 1, 2, 4} (c0 = xi * 0.5 sqrt(N - 1) / ||W_K||_F), held-out seeds 100-102, 12,000 steps, two batches,
       scored against the better of two native SA runs of 64,000 sweeps (as for K4T in round 7).
 CPP   constraint-preserving bSB, run pools for matched wall-clock time (seeds a-b, K = 3 and 4): 3,000, 12,000 and 48,000 steps, eight
       replicas, best plan of the batch, pools sized as in round7 HEATP.
CPTB  the default bSB on the CPT tuning reference.
Usage: python3 bench_round9.py REP|OH2|POLT|CPT|CPTB|CPP a-b  -> results/round9_<part>.json"""
import sys, os, json, math, time, numpy as np, scipy.sparse as sp
from exp_common import *
from refs import build_ref

part = sys.argv[1]; out = {}


def sb_traj(a, J, steps, R, rng, every=50):
    """bSB as sb_batch (dt 0.5) that also returns, per replica, the lowest-energy sign state visited (sampled every `every` steps)."""
    N = len(a); Jm = J * 0.25; h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
    c0 = 0.5 * math.sqrt(N - 1) / (math.sqrt((Jm.multiply(Jm)).sum()) + 1e-12)
    X = rng.uniform(-0.1, 0.1, (N, R)); Y = rng.uniform(-0.1, 0.1, (N, R))
    bestE = np.full(R, np.inf); bestS = np.zeros((N, R)); visited = []
    for t in range(steps):
        at = t / steps
        Y += 0.5 * (-(1.0 - at) * X - c0 * (Jm @ X + h[:, None])); X += 0.5 * Y
        m = np.abs(X) > 1; X[m] = np.sign(X[m]); Y[m] = 0.0
        if (t + 1) % every == 0 or t == steps - 1:
            S = np.where(X >= 0, 1.0, -1.0); E = S.T @ h + 0.5 * np.einsum("ir,ir->r", S, Jm @ S)
            k = E < bestE; bestE[k] = E[k]; bestS[:, k] = S[:, k]
            if t >= steps // 2: visited.append(((S + 1) / 2).astype(np.int64).T.copy())
    return (X > 0).astype(np.int64).T, ((bestS + 1) / 2).astype(np.int64).T, visited


if part == "REP":
    REF = build_ref(); warm(3); warm(4); rows = []
    for K in (3, 4):
        for seed in range(5):
            W, csr, a, J, lam, jj = setup(200, seed, K); rng = np.random.default_rng(11000 + 10 * seed + K); rec = dict(K=K, seed=seed, ref=REF[(K, 200, seed)], runs={})
            for steps in (3000, 12000, 48000):
                b = []
                for _ in range(4):
                    t0 = time.perf_counter(); Xf, Xv, vis = sb_traj(a, J, steps, 8, rng); tt = time.perf_counter() - t0
                    ef = [qubo_energy(x, a, *jj) for x in Xf]; ev = [qubo_energy(x, a, *jj) for x in Xv]
                    ff = [decode_eval(x, K, csr)[0] for x in Xf]; fv = [decode_eval(x, K, csr)[0] for x in Xv]
                    fvis = min(decode_eval(x, K, csr)[0] for V in vis for x in V)
                    b.append(dict(default=float(ff[int(np.argmin(ef))]), all8=float(min(ff)), visited_E=float(fv[int(np.argmin(ev))]),
                                  visited_best=float(min(min(fv), fvis)), same_state=bool(np.array_equal(Xf[int(np.argmin(ef))], Xv[int(np.argmin(ev))])), t=tt))
                rec["runs"][str(steps)] = b
            rows.append(rec); print("REP", K, seed, {s: {k: round((np.median([x[k] for x in v]) / rec["ref"] - 1) * 100, 2) for k in ("default", "all8", "visited_E", "visited_best")} for s, v in rec["runs"].items()}, flush=True)
    out["rows"] = rows

if part == "REP4":     # K = 4: what the best repaired visited state is (raw feasibility, step), with the row numbering and with permuted labels
    REF = build_ref(); warm(4); rows = []
    for seed in range(5):
        W0 = make_instance(200, seed, 4)
        for lab in ("row", "perm"):
            perm = np.arange(200) if lab == "row" else np.random.default_rng(300 + seed).permutation(200)
            W = W0[np.ix_(perm, perm)]; csr = csr_from_W(W); a, J, lam = build_qubo(W, 4, CFG["lam_scale"]); rng = np.random.default_rng(16000 + seed)
            jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64)); res = []
            for b in range(2):
                Xf, Xv, vis = sb_traj(a, J, 12000, 8, rng)
                ef = [qubo_energy(x, a, *jj) for x in Xf]; fdef, phidef = decode_eval(Xf[int(np.argmin(ef))], 4, csr)[:2]
                best = None
                for k, V in enumerate(vis):
                    for x in V:
                        f, phi = decode_eval(x, 4, csr)[:2]
                        if best is None or f < best[0]: best = (float(f), float(phi), 12000 // 2 + 50 * k)
                res.append(dict(default=float(fdef), phi_default=float(phidef), vis_f=best[0], vis_phi=best[1], vis_step=best[2]))
            g = lambda v: round((v / REF[(4, 200, seed)] - 1) * 100, 2)
            rows.append(dict(seed=seed, labels=lab, ref=REF[(4, 200, seed)], runs=res)); print("REP4", seed, lab, [(g(r["default"]), g(r["vis_f"]), round(r["vis_phi"], 2), r["vis_step"]) for r in res], flush=True)
    out["rows"] = rows

if part == "OH2":
    K2 = {r["seed"]: r for r in json.load(open("results/round6_K2.json"))["rows"] if r["sigma"] == 0.15}; warm(3); rows = []
    for seed in range(5):
        W = make_instance(200, seed, 2); csr = csr_from_W(W); a, J, lam = build_qubo(W, 2, CFG["lam_scale"]); jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64))
        rng = np.random.default_rng(12000 + seed); rec = dict(seed=seed, best_round6=K2[seed]["best"], runs={})
        for steps in (3000, 12000, 48000):
            fs = []; ts = []
            for _ in range(4):
                t0 = time.perf_counter(); Xb = sb_batch(a, J, steps, 8, rng, False); es = [qubo_energy(x, a, *jj) for x in Xb]
                f, phi = decode_eval(Xb[int(np.argmin(es))], 2, csr)[:2]; fs.append(dict(f=float(f), phi=float(phi))); ts.append(time.perf_counter() - t0)
            rec["runs"][str(steps)] = dict(rows=fs, t=float(np.mean(ts)))
        rows.append(rec); print("OH2", seed, {s: (round((np.median([x["f"] for x in v["rows"]]) / rec["best_round6"] - 1) * 100, 2), min(x["phi"] for x in v["rows"])) for s, v in rec["runs"].items()}, flush=True)
    out["rows"] = rows; out["lam_scale"] = CFG["lam_scale"]

if part == "POLT":
    import random
    warm(3); rows = []

    def accept_rate(ch0, csr, T0, T1, sweeps, which, seed):
        ip, ix, dt_ = csr; random.seed(seed); ch = ch0.copy(); n = len(ch); acc = 0; tot = 0
        rng_s = range(sweeps) if which == "first" else range(sweeps)
        for s in rng_s:
            T = T0 * (T1 / T0) ** (s / max(sweeps - 1, 1))
            for _ in range(n):
                i = random.randrange(n); c2 = random.randrange(2)
                if c2 >= ch[i]: c2 += 1
                d = 0.0
                for k in range(ip[i], ip[i + 1]):
                    cj = ch[ix[k]]
                    if cj == c2: d += dt_[k]
                    elif cj == ch[i]: d -= dt_[k]
                tot += 1
                if d <= 0 or random.random() < math.exp(-d / T): ch[i] = c2; acc += 1
        return acc / tot

    for seed in range(5):
        W, csr, a, J, lam, jj = setup(200, seed, 3); rng = np.random.default_rng(13000 + seed)
        Xb = sb_batch(a, J, 3000, 8, rng, False); es = [qubo_energy(x, a, *jj) for x in Xb]; ch_b, _ = repair(Xb[int(np.argmin(es))].astype(np.int64), 3, *csr)
        ch_r = rng.integers(0, 3, 200).astype(np.int64); rec = dict(seed=seed)
        for name, ch0, T0 in (("frozen_from_random", ch_r, 2.0), ("polish_from_bsb", ch_b, 0.3)):
            ts = []
            for k in range(3):
                t0 = time.perf_counter(); sa_native(ch0.copy(), 3, *csr, T0, 0.02, 64000, 500 + k); ts.append(time.perf_counter() - t0)
            # acceptance over the whole schedule, estimated on 2,000 sweeps of the same geometric schedule (first and last 1,000 by time)
            rec[name] = dict(t=float(np.median(ts)), acc=accept_rate(ch0, csr, T0, 0.02, 2000, "all", 7))
        rows.append(rec); print("POLT", seed, rec, flush=True)
    out["rows"] = rows


def cp_qubo_parts(W, K):
    n = W.shape[0]; Wk = sp.kron(sp.csr_matrix(W), sp.identity(K)).tocsr()     # x index i*K + c
    return n, Wk


def proj_simplex(V):
    """Euclidean projection of each row of V (shape m x K) on the probability simplex."""
    U = -np.sort(-V, axis=1); css = np.cumsum(U, axis=1) - 1.0; k = np.arange(1, V.shape[1] + 1)
    cond = U - css / k > 0; r = cond.sum(1); theta = css[np.arange(len(V)), r - 1] / r
    return np.maximum(V - theta[:, None], 0.0)


def cp_bsb(Wk, n, K, steps, R, rng, dt, xi):
    N = n * K; c0 = xi * 0.5 * math.sqrt(N - 1) / math.sqrt(Wk.multiply(Wk).sum())
    X = proj_simplex((1.0 / K + rng.uniform(-0.05, 0.05, (n * R, K))))                  # rows: (replica r, cell i)
    X = X.reshape(R, n, K).transpose(1, 2, 0).reshape(N, R)
    Y = rng.uniform(-0.05, 0.05, (N, R)); Y = Y - Y.reshape(n, K, R).mean(1, keepdims=True).repeat(K, 1).reshape(N, R)
    for t in range(steps):
        at = t / steps
        F = -(1.0 - at) * (X - 1.0 / K) - c0 * (Wk @ X)
        F = (F.reshape(n, K, R) - F.reshape(n, K, R).mean(1, keepdims=True)).reshape(N, R)
        Y += dt * F; X += dt * Y
        Z = X.reshape(n, K, R).transpose(2, 0, 1).reshape(R * n, K)
        Zp = proj_simplex(Z); hit = Zp <= 0.0
        X = Zp.reshape(R, n, K).transpose(1, 2, 0).reshape(N, R)
        Yr = Y.reshape(n, K, R).transpose(2, 0, 1).reshape(R * n, K); Yr[hit] = 0.0
        Y = Yr.reshape(R, n, K).transpose(1, 2, 0).reshape(N, R)
    return X.reshape(n, K, R).argmax(1).T                                                # R x n channels


if part == "CPT":
    warm(3); warm(4); res = {}
    for K in (3, 4):
        sa_ref = {}
        for s in (100, 101, 102):
            W = make_instance(200, s, K); csr = csr_from_W(W)
            sa_ref[s] = min(sa_native(np.random.default_rng(s + k).integers(0, K, 200), K, *csr, 2.0, 0.02, 64000, 900 + k)[1] for k in range(2))
        for dt in (0.25, 0.5, 1.0):
            for xi in (0.5, 1.0, 2.0, 4.0):
                v = []
                for s in (100, 101, 102):
                    W = make_instance(200, s, K); csr = csr_from_W(W); n, Wk = cp_qubo_parts(W, K); rng = np.random.default_rng(14000 + s)
                    for b in range(2):
                        C = cp_bsb(Wk, n, K, 12000, 8, rng, dt, xi); f = min(float(interference(c.astype(np.int64), *csr)) for c in C)
                        v.append((f / sa_ref[s] - 1) * 100)
                key = f"K{K}_dt{dt}_xi{xi}"; res[key] = dict(K=K, dt=dt, xi=xi, gaps=v, mean_gap=float(np.mean(v))); print("CPT", key, round(float(np.mean(v)), 2), flush=True)
    out = res

if part == "CPTB":      # the default bSB scored on the same held-out reference as CPT (same SA seeds, 12,000 steps, two batches)
    warm(3); warm(4); res = {}
    for K in (3, 4):
        v = []
        for s in (100, 101, 102):
            W = make_instance(200, s, K); csr = csr_from_W(W)
            ref = min(sa_native(np.random.default_rng(s + k).integers(0, K, 200), K, *csr, 2.0, 0.02, 64000, 900 + k)[1] for k in range(2))
            a, J, lam = build_qubo(W, K, CFG["lam_scale"]); jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64)); rng = np.random.default_rng(14000 + s)
            for b in range(2):
                Xb = sb_batch(a, J, 12000, 8, rng, False); es = [qubo_energy(x, a, *jj) for x in Xb]
                v.append((decode_eval(Xb[int(np.argmin(es))], K, csr)[0] / ref - 1) * 100)
        res[f"K{K}"] = dict(gaps=v, mean_gap=float(np.mean(v))); print("CPTB", K, round(float(np.mean(v)), 2), flush=True)
    out = res

if part == "CPP":
    G = json.load(open("results/round9_CPT.json")); a_, b_ = map(int, sys.argv[2].split("-")); warm(3); warm(4); rows = []
    KS = [int(sys.argv[3])] if len(sys.argv) > 3 else [3, 4]; jl = f"results/round9_CPP_K{''.join(map(str, KS))}.jsonl"
    done_ = {(json.loads(l)["K"], json.loads(l)["seed"]) for l in open(jl)} if os.path.exists(jl) else set()
    best = {K: min((v for v in G.values() if v["K"] == K), key=lambda v: v["mean_gap"]) for K in (3, 4)}
    for K in KS:
        for seed in range(a_, b_):
            if (K, seed) in done_: continue
            W = make_instance(200, seed, K); csr = csr_from_W(W); n, Wk = cp_qubo_parts(W, K); rng = np.random.default_rng(15000 + 10 * seed + K)
            rec = dict(K=K, seed=seed, cfg=dict(dt=best[K]["dt"], xi=best[K]["xi"]), runs={})
            for steps, Tmax in ((3000, 30), (12000, 100), (48000, 100)):
                pool = []; need = None
                while need is None or len(pool) < need:
                    t0 = time.perf_counter(); C = cp_bsb(Wk, n, K, steps, 8, rng, best[K]["dt"], best[K]["xi"])
                    f = min(float(interference(c.astype(np.int64), *csr)) for c in C); pool.append(dict(f=f, t=time.perf_counter() - t0))
                    if need is None: need = max(4, min(60, 2 * int(Tmax // pool[0]["t"])))
                rec["runs"][str(steps)] = pool
            rows.append(rec); append_jsonl(jl, rec); print("CPP", K, seed, {k: len(v) for k, v in rec["runs"].items()}, flush=True)
    out = dict(cfg={K: dict(dt=best[K]["dt"], xi=best[K]["xi"], tune_gap=best[K]["mean_gap"]) for K in (3, 4)}); part = "CPP_cfg"

json.dump(out, open(f"results/round9_{part}.json", "w"), indent=1)

"""Two cheap checks (K = 3, 200 cells, regular layout, seeds 0-4 unless stated).
 C0   coupling scale. The default c0 = 0.5 sqrt(N - 1) / ||J||_F (Eq. (20)) is reported per instance, with the share of ||J||_F^2 that comes
      from the penalty (intra-cell) couplings and the largest eigenvalue of the Ising coupling matrix. bSB (3,000 and 12,000 steps, four
      batches of eight) is then run with a spectral scale c0 = xi / lambda_max(J), xi in {0.5, 1, 2}, against the default.
 POL  polish from a random plan: native SA with the polish schedule (T0 = 0.3 -> 0.02, 1,000 and 4,000 sweeps, three runs) started from the
      decoded bSB plan and from a random plan, on the held-out seeds 100 and 101 at 200 and 500 cells (the pilot of protocol P6).
 K4T  bSB re-tuned at K = 4 with the grid of protocol P5 (step size {0.25, 0.5, 1} x pump exponent {0.5, 1, 2} x coupling multiplier
      {0.5, 1, 2} x penalty scale {0.12, 0.25, 0.5}; 12,000 steps, two batches on held-out seeds 100-102, 200 cells), scored against the
      best of two native SA runs of 64,000 sweeps; the best configuration is then evaluated on seeds 0-4 at 12,000 and 48,000 steps
      (four batches) next to the default. Run as K4T a|b (halves of the grid) and then K4T eval.
 HEAT thermal-fluctuation (heated) bSB in the form of [8]: the momentum update adds heat * y; heat in {0, 0.01, 0.03, 0.06, 0.1, 0.2}
      chosen on held-out seeds 100-102 (K = 3, 200 cells, 12,000 steps, two batches), then evaluated on seeds 0-4 at 3,000, 12,000
      and 48,000 steps (four batches) next to the default.
 HEATP fresh run pools of heated bSB (heat chosen in HEAT) at 12,000, 48,000 and 192,000 steps on seeds 0-4 (K = 3, 200 cells), sized as in
      bench_fresh67.py (at least 2 floor(Tmax / t) batches, Tmax = 30, 100, 100 s), for the matched wall-clock comparison with the pools of
      Table 7. Run as HEATP 0-3 and HEATP 3-5 on the two cores.
Usage: python3 bench_round7.py C0|POL|K4T a|K4T b|K4T eval|HEAT|HEATP a-b -> results/round7_<part>.json"""
import sys, json, math, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
from exp_common import *
from refs import build_ref

part = sys.argv[1]; out = {}


def sb_c0(a, J, steps, R, rng, c0):
    N = len(a); Jm = J * 0.25; h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
    X = rng.uniform(-0.1, 0.1, (N, R)); Y = rng.uniform(-0.1, 0.1, (N, R))
    for t in range(steps):
        at = t / steps; Y += 0.5 * (-(1.0 - at) * X - c0 * (Jm @ X + h[:, None])); X += 0.5 * Y
        m = np.abs(X) > 1; X[m] = np.sign(X[m]); Y[m] = 0.0
    return (X > 0).astype(np.int64).T


if part == "C0":
    REF = build_ref(); warm(3); rows = []
    for seed in range(5):
        W, csr, a, J, lam, jj = setup(200, seed, 3); N = len(a); Jm = J * 0.25
        fro = math.sqrt(Jm.multiply(Jm).sum()); c0 = 0.5 * math.sqrt(N - 1) / fro
        Jc = Jm.tocoo(); intra = (Jc.row // 3) == (Jc.col // 3); pen_share = float((Jc.data[intra] ** 2).sum() / fro ** 2)
        lmax = float(spla.eigsh(Jm.tocsc().astype(float), k=1, which="LA", return_eigenvectors=False)[0])
        rec = dict(seed=seed, c0=c0, fro=fro, pen_share=pen_share, lmax=lmax, ref=REF[(3, 200, seed)], runs={})
        rng = np.random.default_rng(7100 + seed)
        for lab, cc in (("default", c0), ("spec0.5", 0.5 / lmax), ("spec1", 1.0 / lmax), ("spec2", 2.0 / lmax)):
            for steps in (3000, 12000):
                fs = []
                for b in range(4):
                    Xb = sb_c0(a, J, steps, 8, rng, cc); es = [qubo_energy(x, a, *jj) for x in Xb]
                    fs.append(float(decode_eval(Xb[int(np.argmin(es))], 3, csr)[0]))
                rec["runs"][f"{lab}|{steps}"] = dict(c0=cc, f=fs)
        rows.append(rec)
        print("C0", seed, round(c0, 4), round(pen_share, 3), round(lmax, 2), {k: round((np.median(v["f"]) / rec["ref"] - 1) * 100, 1) for k, v in rec["runs"].items()}, flush=True)
    out["rows"] = rows
if part == "POL":
    warm(3); rows = []
    for n in (200, 500):
        for seed in (100, 101):
            W, csr, a, J, lam, jj = setup(n, seed, 3); rng = np.random.default_rng(5)
            Xb = sb_batch(a, J, 4 * 3000, 8, rng, False); es = [qubo_energy(x, a, *jj) for x in Xb]; rr = repair_all(Xb, 3, csr)
            ch0, f0 = rr[int(np.argmin(es))]
            rec = dict(n=n, seed=seed, f_bsb=float(f0))
            for sw in (1000, 4000):
                rec[f"from_bsb_{sw}"] = [float(sa_native(ch0, 3, *csr, 0.3, 0.02, sw, s)[1]) for s in range(3)]
                rec[f"from_random_{sw}"] = [float(sa_native(rng.integers(0, 3, n), 3, *csr, 0.3, 0.02, sw, s)[1]) for s in range(3)]
            rows.append(rec); print("POL", n, seed, round(f0, 2), {k: round(float(np.mean(v)), 2) for k, v in rec.items() if k.startswith("from")}, flush=True)
    out["rows"] = rows
if part == "K4T":
    import itertools
    from bench_round5 import sb2
    sub = sys.argv[2]; warm(4); K = 4
    grid = [dict(dt=dt, gamma=g, c0m=c, rho=r) for dt, g, c, r in itertools.product((0.25, 0.5, 1.0), (0.5, 1.0, 2.0), (0.5, 1.0, 2.0), (0.12, 0.25, 0.5))]
    if sub in ("a", "b"):
        grid = grid[:41] if sub == "a" else grid[41:]
        inst = {}
        for s_ in (100, 101, 102):
            W = make_instance(200, s_, K); csr = csr_from_W(W)
            sa_ = min(float(sa_native(np.random.default_rng(s_).integers(0, K, 200), K, *csr, 0.15, 0.02, 64000, 500 + k)[1]) for k in range(2))
            inst[s_] = (W, csr, sa_)
        res = {}
        for cfg in grid:
            key = "dt%g_g%g_c%g_rho%g" % (cfg["dt"], cfg["gamma"], cfg["c0m"], cfg["rho"]); v = []
            for s_, (W, csr, sa_) in inst.items():
                a, J, lam = build_qubo(W, K, cfg["rho"]); jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64))
                rng = np.random.default_rng(900 + s_)
                for b in range(2):
                    Xb = sb2(a, J, 12000, 8, rng, cfg["dt"], cfg["gamma"], cfg["c0m"]); es = [qubo_energy(x, a, *jj) for x in Xb]
                    v.append(dict(seed=s_, f=float(decode_eval(Xb[int(np.argmin(es))], K, csr)[0]), sa=sa_))
            res[key] = dict(cfg=cfg, runs=v, mean_gap=float(np.mean([(x["f"] / x["sa"] - 1) * 100 for x in v])))
            print("K4T", key, round(res[key]["mean_gap"], 2), flush=True)
        out = res; part = f"K4T_{sub}"
    else:
        G = {**json.load(open("results/round7_K4T_a.json")), **json.load(open("results/round7_K4T_b.json"))}
        # a configuration whose dynamics diverge returns the same clipped state in every batch, and its decoded plan is the decoder-only
        # repair of that state; such configurations are excluded by requiring a raw feasibility of at least 0.95 on held-out seed 100
        W0 = make_instance(200, 100, K); csr0 = csr_from_W(W0); checked = {}
        for g_, k_ in sorted((G[k]["mean_gap"], k) for k in G):
            c = G[k_]["cfg"]; a, J, lam = build_qubo(W0, K, c["rho"]); jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64))
            Xb = sb2(a, J, 12000, 8, np.random.default_rng(1), c["dt"], c["gamma"], c["c0m"]); es = [qubo_energy(x, a, *jj) for x in Xb]
            phi = float(decode_eval(Xb[int(np.argmin(es))], K, csr0)[1]); checked[k_] = phi; print("check", k_, round(g_, 2), "phi", phi, flush=True)
            if phi >= 0.95: best = k_; break
        cfg = G[best]["cfg"]; REF = build_ref(); rows = []
        for seed in range(5):
            W = make_instance(200, seed, K); csr = csr_from_W(W); rec = dict(seed=seed, ref=REF[(4, 200, seed)], runs={})
            for lab, c in (("default", dict(dt=0.5, gamma=1.0, c0m=1.0, rho=0.25)), ("tuned", cfg)):
                a, J, lam = build_qubo(W, K, c["rho"]); jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64))
                rng = np.random.default_rng(950 + seed)
                for steps in (12000, 48000):
                    fs = []
                    for b in range(4):
                        Xb = sb2(a, J, steps, 8, rng, c["dt"], c["gamma"], c["c0m"]); es = [qubo_energy(x, a, *jj) for x in Xb]
                        fs.append(float(decode_eval(Xb[int(np.argmin(es))], K, csr)[0]))
                    rec["runs"][f"{lab}|{steps}"] = fs
            rows.append(rec); print("K4T eval", seed, {k: round((np.median(v) / rec["ref"] - 1) * 100, 1) for k, v in rec["runs"].items()}, flush=True)
        out = dict(best=best, cfg=cfg, phi_checked=checked, grid_best=sorted(((G[k]["mean_gap"], k) for k in G))[:5], default_grid=G.get("dt0.5_g1_c1_rho0.25", {}).get("mean_gap"), rows=rows)
        part = "K4T_eval"
if part == "HEAT":
    exec(open("bench_round3.py").read().split("def decode_best")[0].split("import exp_common")[0])
    src = open("bench_round3.py").read(); i0 = src.index("def sb_custom"); i1 = src.index("def decode_best"); exec(src[i0:i1])
    warm(3); REF = build_ref(); K = 3; tune = {}
    inst = {}
    for s_ in (100, 101, 102):
        W, csr, a, J, lam, jj = setup(200, s_, K)
        sa_ = min(float(sa_native(np.random.default_rng(s_).integers(0, K, 200), K, *csr, 0.6, 0.02, 64000, 600 + k)[1]) for k in range(2))
        inst[s_] = (W, csr, a, J, jj, sa_)
    for heat in (0.0, 0.01, 0.03, 0.06, 0.1, 0.2):
        v = []
        for s_, (W, csr, a, J, jj, sa_) in inst.items():
            rng = np.random.default_rng(980 + s_)
            for b in range(2):
                Xb = sb_custom(a, J, 12000, 8, rng, heat=heat); es = [qubo_energy(x, a, *jj) for x in Xb]
                v.append((float(decode_eval(Xb[int(np.argmin(es))], K, csr)[0]) / sa_ - 1) * 100)
        tune[str(heat)] = float(np.mean(v)); print("HEAT tune", heat, round(tune[str(heat)], 2), flush=True)
    hb = float(min(tune, key=tune.get)); rows = []
    for seed in range(5):
        W, csr, a, J, lam, jj = setup(200, seed, K); rec = dict(seed=seed, ref=REF[(3, 200, seed)], runs={})
        for lab, ht in (("default", 0.0), ("heated", hb)):
            rng = np.random.default_rng(990 + seed)
            for steps in (3000, 12000, 48000):
                fs = []
                for b in range(4):
                    Xb = sb_custom(a, J, steps, 8, rng, heat=ht); es = [qubo_energy(x, a, *jj) for x in Xb]
                    fs.append(float(decode_eval(Xb[int(np.argmin(es))], K, csr)[0]))
                rec["runs"][f"{lab}|{steps}"] = fs
        rows.append(rec); print("HEAT eval", seed, {k: round((np.median(v) / rec["ref"] - 1) * 100, 1) for k, v in rec["runs"].items()}, flush=True)
    out = dict(tune=tune, heat=hb, rows=rows)
if part == "HEATP":
    src = open("bench_round3.py").read(); i0 = src.index("def sb_custom"); i1 = src.index("def decode_best"); exec(src[i0:i1])
    hb = json.load(open("results/round7_HEAT.json"))["heat"]; a_, b_ = map(int, sys.argv[2].split("-")); warm(3); rows = []
    for seed in range(a_, b_):
        W, csr, a, J, lam, jj = setup(200, seed, 3); rng = np.random.default_rng(9900 + seed); rec = dict(seed=seed, heat=hb, runs={})
        for steps, Tmax in ((12000, 30), (48000, 100), (192000, 100)):
            pool = []; need = None
            while need is None or len(pool) < need:
                t0 = time.perf_counter(); Xb = sb_custom(a, J, steps, 8, rng, heat=hb); es = [qubo_energy(x, a, *jj) for x in Xb]
                f, phi = decode_eval(Xb[int(np.argmin(es))], 3, csr)[:2]; pool.append(dict(f=float(f), phi=float(phi), t=time.perf_counter() - t0))
                if need is None: need = max(4, 2 * int(Tmax // pool[0]["t"]))
            rec["runs"][str(steps)] = pool
        rows.append(rec); print("HEATP", seed, {k: len(v) for k, v in rec["runs"].items()}, flush=True)
    out = dict(rows=rows); part = f"HEATP_{a_}_{b_}"
json.dump(out, open(f"results/round7_{part}.json", "w"), indent=1)

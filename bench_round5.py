"""Round-5 checks: (A) common-seed penalty sweep at 200 cells (K=3, 12 s per solver and instance, restart counts);
(B) cross-validation of bSB tuning: top-5 stage-1 configurations + default re-scored on 10 fresh held-out instances
(n=200, K=3, seeds 110-119; 12k steps, 2 batches each), gaps vs per-instance best of long native SA and all batches;
(C) sensitivity to path-loss exponent eta and interference radius d_max: reuse, native SA (64k sweeps, 3 runs) and
reuse+SA on 5 instances (K=3, n=200) for eta in {3, 3.5, 4} x d_max in {2.0, 2.6, 3.5}."""
import json, math, time, itertools, numpy as np
from exp_common import *
from bench_ran_qubo import sa_qubo
from bench_reuse_polish import best_reuse
import importlib.util
spec = importlib.util.spec_from_file_location("tb", "/root/.claude/uploads/7e3b9e47-e110-5520-9822-bb0ac6b242f3/7914bcc9-tune_bsb.py")


def sb2(a, J, steps, R, rng, dt, gamma, c0m):
    N = len(a); Jm = J * 0.25; h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
    c0 = c0m * 0.5 * math.sqrt(N - 1) / (math.sqrt(Jm.multiply(Jm).sum()) + 1e-12)
    X = rng.uniform(-0.1, 0.1, (N, R)); Y = rng.uniform(-0.1, 0.1, (N, R))
    for t in range(steps):
        at = (t / steps) ** gamma; Y += dt * (-(1.0 - at) * X - c0 * (Jm @ X + h[:, None])); X += dt * Y
        m = np.abs(X) > 1; X[m] = np.sign(X[m]); Y[m] = 0.0
    return (X > 0).astype(np.int64).T


def part_a():
    out = "results/pen200.jsonl"; done = done_keys(out); K = 3; n = 200; B = 12.0
    for seed in range(5):
        if (n, seed, K) in done: continue
        W = make_instance(n, seed, K); csr = csr_from_W(W); N = n * K; res = {}
        def loop(fn):
            rng = np.random.default_rng(8000 + seed); t0 = time.perf_counter(); best = 1e300; k = 0; fe = []
            while time.perf_counter() - t0 < B:
                e, f = fn(rng, k); best = min(best, e); fe.append(f); k += 1
            return dict(E=float(best), restarts=k, feas=float(np.mean(fe)))
        res["SA-native"] = loop(lambda rng, k: (sa_native(rng.integers(0, K, n), K, *csr, 2.0, 0.02, CFG["native_sweeps"], 80000 + k)[1], 1.0))
        for ls in (1.0, 0.5, 0.25, 0.12, 0.06):
            a, J, lam = build_qubo(W, K, ls); jp, jx, jd = J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64)
            res[f"SA-QUBO|{ls}"] = loop(lambda rng, k: decode_eval(sa_qubo(rng.integers(0, 2, N), a, jp, jx, jd, 0.5 * lam, 0.15, CFG["sweeps"], 81000 + k)[0], K, csr))
            def sbf(rng, k):
                Xb = sb_batch(a, J, CFG["sb_steps"], CFG["sb_R"], rng, False); es = [qubo_energy(z, a, jp, jx, jd) for z in Xb]
                return decode_eval(Xb[int(np.argmin(es))], K, csr)
            res[f"bSB|{ls}"] = loop(sbf)
        append_jsonl(out, dict(n=n, seed=seed, K=K, budget=B, res=res)); print("A", seed, flush=True)


def part_b():
    g = json.load(open("paper_data/tuning_grid.json"))["grid"]
    score = {k: np.mean([r[1] for r in v["runs"]]) for k, v in g.items()}
    top = sorted(score, key=score.get)[:5]
    cfgs = {k: g[k]["cfg"] for k in top}; cfgs["default"] = dict(dt=0.5, gamma=1.0, c0m=1.0, rho=0.25)
    out = "results/bsb_cv.jsonl"; K = 3; n = 200
    for seed in range(110, 120):
        W = make_instance(n, seed, K); csr = csr_from_W(W); rng = np.random.default_rng(4200 + seed)
        best = min(sa_native(rng.integers(0, K, n), K, *csr, 2.0, 0.02, 64000, 600 + k)[1] for k in range(2)); res = {}
        for name, c in cfgs.items():
            a, J, lam = build_qubo(W, K, c["rho"]); jp, jx, jd = J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64)
            fs = []
            for b in range(2):
                Xb = sb2(a, J, 12000, 8, rng, c["dt"], c["gamma"], c["c0m"]); es = [qubo_energy(x, a, jp, jx, jd) for x in Xb]
                fs.append(float(decode_eval(Xb[int(np.argmin(es))], K, csr)[0]))
            res[name] = fs
        append_jsonl(out, dict(n=n, seed=seed, K=K, sa_best=float(best), tuning_score={k: float(score[k]) for k in top}, res=res)); print("B", seed, flush=True)


def make_eta(n, seed, K, eta, dmax):
    rng = np.random.default_rng(seed); cols = int(math.ceil(math.sqrt(n * 1.15))); pts = []; r = 0
    while len(pts) < n:
        for c in range(cols): pts.append((c + 0.5 * (r % 2), r * math.sqrt(3) / 2))
        r += 1
    pts = np.array(pts[:n]) + rng.normal(0, 0.15, (n, 2))
    d = np.linalg.norm(pts[:, None] - pts[None], axis=2); np.fill_diagonal(d, np.inf)
    W = np.where(d <= dmax, np.maximum(d, 0.5) ** -eta, 0.0); W = np.triu(W, 1); return W + W.T


def part_c():
    out = "results/eta_dmax.jsonl"; K = 3; n = 200
    for eta, dmax in itertools.product((3.0, 3.5, 4.0), (2.0, 2.6, 3.5)):
        for seed in range(5):
            W = make_eta(n, seed, K, eta, dmax); csr = csr_from_W(W); rng = np.random.default_rng(300 + seed)
            p = best_reuse(n, K, csr); f_re = float(interference(p, *csr))
            sa = [float(sa_native(rng.integers(0, K, n), K, *csr, 2.0, 0.02, 64000, 700 + k)[1]) for k in range(3)]
            rs = float(sa_native(p, K, *csr, 0.3, 0.02, 64000, 900)[1])
            append_jsonl(out, dict(eta=eta, dmax=dmax, seed=seed, reuse=f_re, sa=sa, reuse_sa=rs))
        print("C", eta, dmax, flush=True)


if __name__ == "__main__":
    warm(3); part_a(); part_b(); part_c(); print("DONE", flush=True)

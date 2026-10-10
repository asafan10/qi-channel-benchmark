"""Part T: tuned single-bit-flip SA on the one-hot QUBO. Schedule grid (T0, T1) chosen per penalty on held-out seeds 100-101 at 200 cells
        (64,000 sweeps, two runs), then evaluated on seeds 0-9 (two runs) and compared, paired, with bSB at matched wall-clock time
        (bSB batches of 3,000 steps and 8 replicas until the SA time is used; best decoded plan).  -> results/qubo_tuned.jsonl
part P: planted 3-colourings at 200 cells: the layout and weights of make_instance with every edge between cells of equal planted
        colour removed, so the optimum is 0. bSB (default coupling scale, and the scale computed from the inter-cell couplings only)
        at 3,000/12,000/48,000 steps, and native SA, all decoded and repaired.  -> results/planted.jsonl"""
import sys, time, math, json
import numpy as np, scipy.sparse as sp
from exp_common import *

GRID = [("own", None, 0.15), ("0.5lam-0.02", None, 0.02), ("0.5-0.02", 0.5, 0.02), ("1.0-0.02", 1.0, 0.02), ("2.0-0.02", 2.0, 0.02)]


def sa_run(W, K, rho, g, S, sd):
    a, J, lam = build_qubo(W, K, rho); jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64))
    T0 = 0.5 * lam if g[1] is None else g[1]
    x0 = np.zeros(len(a), dtype=np.int64); ch = np.random.default_rng(sd).integers(0, K, W.shape[0]); x0[np.arange(W.shape[0]) * K + ch] = 1
    t0 = time.perf_counter(); bx, _ = sa_qubo(x0, a, *jj, T0, g[2], S, sd); return bx, time.perf_counter() - t0, (a, J, jj)


def part_t():
    K, n = 3, 200; out = "results/qubo_tuned.jsonl"; best = {}
    for rho in (0.25, 0.12, 0.06):
        score = {}
        for name, *g in GRID:
            cs = []
            for seed in (100, 101):
                W = make_instance(n, seed, K); csr = csr_from_W(W)
                for r in range(2): cs.append(decode_eval(sa_run(W, K, rho, (name, *g), 64000, 5000 + 10 * seed + r)[0], K, csr)[0])
            score[name] = float(np.mean(cs))
        best[rho] = min(score, key=score.get); print("T select", rho, score, best[rho], flush=True)
        g = next(x for x in GRID if x[0] == best[rho])
        for seed in range(10):
            W = make_instance(n, seed, K); csr = csr_from_W(W); fs = []; ts = 0.0
            for r in range(2):
                bx, t, (a, J, jj) = sa_run(W, K, rho, g, 64000, 7000 + 10 * seed + r); fs.append(float(decode_eval(bx, K, csr)[0])); ts += t
            rng = np.random.default_rng(8000 + seed); fb = []; tb = 0.0
            while tb < ts:
                t0 = time.perf_counter(); Xb = sb_batch(a, J, 3000, 8, rng, False); es = [qubo_energy(x, a, *jj) for x in Xb]
                fb.append(float(decode_eval(Xb[int(np.argmin(es))], K, csr)[0])); tb += time.perf_counter() - t0
            append_jsonl(out, dict(n=n, seed=seed, K=K, rho=rho, sched=best[rho], scores=score, sa=fs, sa_t=ts, bsb=fb, bsb_t=tb))
        print("T done", rho, flush=True)


def sb_c0(a, J, steps, R, rng, c0mode):
    N = len(a); Jm = J * 0.25; h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
    if c0mode == "inter":
        K = 3; C = J.tocoo(); keep = (C.row // K) != (C.col // K); Ji = sp.coo_matrix((C.data[keep] * 0.25, (C.row[keep], C.col[keep])), shape=J.shape)
        scale = math.sqrt((Ji.multiply(Ji)).sum()) + 1e-12
    else:
        scale = math.sqrt((Jm.multiply(Jm)).sum()) + 1e-12
    c0 = 0.5 * math.sqrt(N - 1) / scale
    X = rng.uniform(-0.1, 0.1, (N, R)); Y = rng.uniform(-0.1, 0.1, (N, R))
    for t in range(steps):
        at = t / steps; force = -(Jm @ X + h[:, None]); Y += 0.5 * (-(1.0 - at) * X + c0 * force); X += 0.5 * Y
        m = np.abs(X) > 1; X[m] = np.sign(X[m]); Y[m] = 0.0
    return (X > 0).astype(np.int64).T


def part_p():
    K, n = 3, 200; out = "results/planted.jsonl"
    for seed in range(5):
        W0 = make_instance(n, seed, K); col = np.random.default_rng(900 + seed).integers(0, K, n)
        W = W0 * (col[:, None] != col[None, :]); csr = csr_from_W(W)
        assert interference(col.astype(np.int64), *csr) == 0.0
        a, J, lam = build_qubo(W, K, 0.25); jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64))
        rec = dict(n=n, seed=seed, K=K, scale=float(W0[np.triu_indices(n, 1)].sum() / 3), res={})
        rng = np.random.default_rng(seed)
        for mode in ("default", "inter"):
            for steps in (3000, 12000, 48000):
                fs = []
                for b in range(4):
                    Xb = sb_c0(a, J, steps, 8, rng, mode); fs += [float(decode_eval(x, K, csr)[0]) for x in Xb]
                rec["res"][f"bSB-{mode}|{steps}"] = fs
        rec["res"]["SA-native|64000"] = [float(sa_native(rng.integers(0, K, n), K, *csr, 2.0, 0.02, 64000, 300 + k)[1]) for k in range(4)]
        rec["res"]["SA-native-tuned|64000"] = [float(sa_native(rng.integers(0, K, n), K, *csr, 0.6, 0.02, 64000, 400 + k)[1]) for k in range(4)]
        append_jsonl(out, rec); print("P", seed, {k: (round(min(v), 3), round(float(np.median(v)), 3)) for k, v in rec["res"].items()}, flush=True)


if __name__ == "__main__":
    warm(3)
    if sys.argv[1] == "T": part_t()
    else: part_p()

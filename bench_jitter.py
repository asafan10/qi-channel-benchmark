"""Layout sweep: where does the periodic reuse pattern stop dominating?
Same generator as make_instance (jittered hexagonal lattice, path-loss weights, d <= 2.6, exponent 3.5) with jitter
sigma in {0.15, 0.3, 0.5}, plus a Poisson layout (uniform random positions over the same area). n = 200, K in {3, 4}, seeds 0-7.
Per instance: best reuse pattern (lattice index; for Poisson, the pattern of the nearest lattice site), native SA from random
(64k and 256k sweeps, 3 runs each), native SA polish from reuse (T0 0.3 and 0.15 -> 0.02, 64k sweeps), bSB at baseline steps
(4 batches) decoded, bSB warm-started from the reuse pattern (positions +-0.5) + decode, and each bSB plan polished (16k sweeps).
-> results/jitter.jsonl"""
import sys, json, math, itertools, numpy as np, scipy.sparse as sp
from exp_common import *
from bench_reuse_polish import best_reuse


def make_layout(n, seed, K, sigma, poisson=False):
    rng = np.random.default_rng(10000 + seed + int(1000 * sigma) + (77777 if poisson else 0))
    cols = int(math.ceil(math.sqrt(n * 1.15))); pts = []; r = 0
    while len(pts) < n:
        for c in range(cols): pts.append((c + 0.5 * (r % 2), r * math.sqrt(3) / 2))
        r += 1
    base = np.array(pts[:n])
    if poisson:
        lo, hi = base.min(0), base.max(0); P = rng.uniform(lo, hi, (n, 2))
        # relabel each random point by its nearest free lattice site so the reuse pattern is defined by location
        order = np.argsort(np.linalg.norm(P[:, None] - base[None], axis=2).min(1)); used = set(); lab = np.empty(n, int)
        for i in order:
            dd = np.linalg.norm(base - P[i], axis=1); dd[list(used)] = np.inf; j = int(np.argmin(dd)); used.add(j); lab[i] = j
        pts = np.empty_like(P); pts[lab] = P
    else:
        pts = base + rng.normal(0, sigma, (n, 2))
    d = np.linalg.norm(pts[:, None] - pts[None], axis=2); np.fill_diagonal(d, np.inf)
    W = np.where(d <= 2.6, np.maximum(d, 0.5) ** -3.5, 0.0); W = np.triu(W, 1); return W + W.T


def bsb_warm(a, J, steps, R, rng, s0):
    """our bSB with positions initialised at 0.5 * reuse spins (+ small noise)"""
    N = len(a); Jm = J * 0.25; h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
    c0 = 0.5 * math.sqrt(N - 1) / (math.sqrt((Jm.multiply(Jm)).sum()) + 1e-12)
    X = 0.5 * s0[:, None] + rng.uniform(-0.05, 0.05, (N, R)); Y = rng.uniform(-0.1, 0.1, (N, R))
    for t in range(steps):
        at = t / steps; Y += 0.5 * (-(1.0 - at) * X - c0 * (Jm @ X + h[:, None])); X += 0.5 * Y
        m = np.abs(X) > 1; X[m] = np.sign(X[m]); Y[m] = 0.0
    return (X > 0).astype(np.int64).T


if __name__ == "__main__":
    out = "results/jitter.jsonl"; done = set()
    import os
    if os.path.exists(out):
        for l in open(out):
            d = json.loads(l); done.add((d["layout"], d["K"], d["seed"]))
    warm(3); warm(4); n = 200
    for layout in ("0.15", "0.3", "0.5", "poisson"):
        for K in (3, 4):
            for seed in range(8):
                if (layout, K, seed) in done: continue
                W = make_layout(n, seed, K, 0.0 if layout == "poisson" else float(layout), layout == "poisson")
                csr = csr_from_W(W); a, J, lam = build_qubo(W, K, CFG["lam_scale"]); jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64))
                rng = np.random.default_rng(555 + seed); res = {}
                p = best_reuse(n, K, csr); res["reuse"] = [float(interference(p, *csr))]
                res["SA x64"] = [float(sa_native(rng.integers(0, K, n), K, *csr, 2.0, 0.02, 64000, 1000 + k)[1]) for k in range(3)]
                res["SA x256"] = [float(sa_native(rng.integers(0, K, n), K, *csr, 2.0, 0.02, 256000, 2000 + k)[1]) for k in range(3)]
                res["reuse+SA"] = [float(sa_native(p, K, *csr, T0, 0.02, 64000, 3000 + k)[1]) for T0 in (0.3, 0.15) for k in range(2)]
                xs = np.zeros(n * K, np.int64); xs[np.arange(n) * K + p] = 1; s0 = 2.0 * xs - 1
                for name, gen in (("bSB", lambda: sb_batch(a, J, CFG["sb_steps"], CFG["sb_R"], rng, False)),
                                  ("bSB from reuse", lambda: bsb_warm(a, J, CFG["sb_steps"], CFG["sb_R"], rng, s0))):
                    fs = []; fp = []
                    for b in range(4):
                        Xb = gen(); es = [qubo_energy(x, a, *jj) for x in Xb]; ch, _ = repair(Xb[int(np.argmin(es))].astype(np.int64), K, *csr)
                        fs.append(float(interference(ch, *csr))); fp.append(float(sa_native(ch, K, *csr, 0.3, 0.02, 16000, 4000 + b)[1]))
                    res[name] = fs; res[name + "+polish"] = fp
                append_jsonl(out, dict(layout=layout, K=K, seed=seed, n=n, res=res))
                best = min(min(v) for v in res.values())
                print(layout, K, seed, {k: round((min(v) / best - 1) * 100, 1) for k, v in res.items()}, flush=True)

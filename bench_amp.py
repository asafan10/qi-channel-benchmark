"""Does replica collapse contribute to the bSB plateau? Initial amplitude of the oscillators.
Our bSB starts positions and momenta in U(-A, A) with A = 0.1. Here A in {0.1, 0.5, 1.0}, everything else unchanged (Eq. (19), dt = 0.5,
linear pump, inelastic walls, 8 replicas per batch, 4 batches).
 G  G22 (Table B6) at 3,000 and 30,000 steps: best and median cut of the 32 replicas, distinct cuts, mean pairwise Hamming distance up to a flip.
 L  irregular layout (jitter 0.3, K = 3, 200 cells, rho = 0.25, seeds 0-4) at 3,000 steps: gap of the best of 32 replicas (best-energy
    replica of each batch, decoded and repaired) to the reference of Table 15, raw feasibility, and replica diversity.
Usage: python3 bench_amp.py G|L -> results/amp_<part>.json"""
import sys, json, math, numpy as np
import exp_common
from exp_common import *
from bench_ran_qubo import qubo_energy
from bench_jitter import make_layout
AMPS = (0.1, 0.5, 1.0)


def sb_amp(a, J, steps, R, rng, A, dt=0.5):
    N = len(a); Jm = J * 0.25; h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
    c0 = 0.5 * math.sqrt(N - 1) / (math.sqrt((Jm.multiply(Jm)).sum()) + 1e-12)
    X = rng.uniform(-A, A, (N, R)); Y = rng.uniform(-A, A, (N, R))
    for t in range(steps):
        at = t / steps; Y += dt * (-(1.0 - at) * X - c0 * (Jm @ X + h[:, None])); X += dt * Y
        m = np.abs(X) > 1; X[m] = np.sign(X[m]); Y[m] = 0.0
    return (X > 0).astype(np.int64).T


def diversity(X):
    S = 2.0 * X - 1.0; D = [(S[i] != S[j]).mean() for i in range(len(S)) for j in range(i + 1, len(S))]
    return float(np.mean([min(d, 1 - d) for d in D]))


part = sys.argv[1]; out = {}
if part == "G":
    from gset_check import read_gset, cut
    Wm = read_gset("gset/G22"); deg = np.asarray(Wm.sum(1)).ravel(); a = -deg; J = (2.0 * Wm).tocsr()
    for A in AMPS:
        for steps in (3000, 30000):
            rng = np.random.default_rng(5); X = np.vstack([sb_amp(a, J, steps, 8, rng, A) for _ in range(4)])
            cuts = [float(cut(Wm, 2.0 * x - 1.0)) for x in X]
            out[f"{A}|{steps}"] = dict(best=max(cuts), median=float(np.median(cuts)), gap_best=(1 - max(cuts) / 13359) * 100,
                                       gap_median=(1 - float(np.median(cuts)) / 13359) * 100, distinct_cuts=len(set(cuts)), mean_hamming=diversity(X))
            print("G", A, steps, out[f"{A}|{steps}"], flush=True)
if part == "L":
    exp_common.make_instance = lambda n, seed, K: make_layout(n, seed, K, 0.3)
    src = open("analyse_sigma_ext.py").read(); exec(src[:src.index("LB = {")])
    for A in AMPS:
        rows = []
        for seed in range(5):
            W, csr, a, J, lam, jj = setup(200, seed, 3); rng = np.random.default_rng(900 + seed); fs, phis, divs = [], [], []
            for b in range(4):
                X = sb_amp(a, J, 3000, 8, rng, A); es = [qubo_energy(x, a, *jj) for x in X]; x = X[int(np.argmin(es))]
                ch, _ = repair(x.astype(np.int64), 3, *csr); fs.append(float(interference(ch, *csr)))
                phis.append(float(np.mean(x.reshape(-1, 3).sum(1) == 1))); divs.append(diversity(X))
            rows.append(dict(seed=seed, gap_best=(min(fs) / REF[seed] - 1) * 100, gap_median_batch=float(np.median([(f / REF[seed] - 1) * 100 for f in fs])),
                             phi=float(np.mean(phis)), hamming=float(np.mean(divs))))
        out[str(A)] = dict(rows=rows, median_gap_best=float(np.median([r["gap_best"] for r in rows])),
                           median_hamming=float(np.median([r["hamming"] for r in rows])), phi=float(np.mean([r["phi"] for r in rows])))
        print("L", A, {k: v for k, v in out[str(A)].items() if k != "rows"}, flush=True)
json.dump(out, open(f"results/amp_{part}.json", "w"), indent=1)

"""Gate check: does our bSB implementation behave as published SB does on unconstrained Max-Cut?
Instances of the G-set families, generated here (G-set itself is not reachable from this machine):
  ER800  : 800 nodes, edge probability 6%, unit weights      (family of G1-G5)
  ER2000 : 2000 nodes, edge probability 1%, unit weights     (family of G22-G26)
  TOR800 : 800-node toroidal 2D grid (20 x 40), +-1 weights  (family of G11-G13)
  G-set files, if present as gset/G<k>, are used too and compared with their best-known cuts.
Our bSB (bench_ran_qubo.sb_batch, unchanged, same c0 rule and pump) runs on the Max-Cut QUBO a_i = -deg_i, J_ij = 2 w_ij,
whose Ising form has zero field. Reference: best of long native single-flip SA runs.
Output: results/gset_check.jsonl"""
import json, os, sys, time, math
import numpy as np, scipy.sparse as sp
from numba import njit
from bench_ran_qubo import sb_batch

GSET_BEST = {1: 11624, 2: 11620, 3: 11622, 11: 564, 12: 556, 22: 13359, 23: 13344, 32: 1410, 43: 6660}


def er(n, p, seed):
    rng = np.random.default_rng(seed); U = sp.triu(sp.random(n, n, density=p, random_state=rng, data_rvs=lambda k: np.ones(k)), 1)
    return (U + U.T).tocsr()


def torus(r, c, seed):
    rng = np.random.default_rng(seed); n = r * c; rows, cols, vals = [], [], []
    for i in range(r):
        for j in range(c):
            u = i * c + j
            for v in (i * c + (j + 1) % c, ((i + 1) % r) * c + j):
                rows.append(min(u, v)); cols.append(max(u, v)); vals.append(rng.choice([-1.0, 1.0]))
    U = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)); return (U + U.T).tocsr()


def read_gset(path):
    with open(path) as fh:
        n, m = map(int, fh.readline().split()); r, c, v = [], [], []
        for _ in range(m):
            i, j, w = fh.readline().split(); r.append(int(i) - 1); c.append(int(j) - 1); v.append(float(w))
    U = sp.coo_matrix((v, (r, c)), shape=(n, n)); return (U + U.T).tocsr()


def cut(Wm, s):            # s in {-1,+1}
    return 0.25 * (Wm.sum() - s @ (Wm @ s))


@njit(cache=True)
def sa_cut(ip, ix, dw, n, sweeps, T0, T1, seed):
    np.random.seed(seed); s = np.where(np.random.random(n) < 0.5, 1.0, -1.0)
    best = -1e18; cur = 0.0
    for i in range(n):
        for k in range(ip[i], ip[i + 1]):
            if ix[k] > i and s[i] != s[ix[k]]: cur += dw[k]
    best = cur
    for t in range(sweeps):
        T = T0 * (T1 / T0) ** (t / max(1, sweeps - 1))
        for _ in range(n):
            i = np.random.randint(n); g = 0.0
            for k in range(ip[i], ip[i + 1]): g += dw[k] * s[i] * s[ix[k]]
            if g >= 0 or np.random.random() < math.exp(g / T):
                s[i] = -s[i]; cur += g
                if cur > best: best = cur
    return best


def run(name, Wm, known=None):
    n = Wm.shape[0]; deg = np.asarray(Wm.sum(1)).ravel(); a = -deg; J = (2.0 * Wm).tocsr()
    ip, ix, dw = Wm.indptr.astype(np.int64), Wm.indices.astype(np.int64), Wm.data.astype(np.float64)
    t0 = time.perf_counter(); ref = max(sa_cut(ip, ix, dw, n, 20000, 3.0, 0.05, 11 + k) for k in range(3)); t_sa = time.perf_counter() - t0
    if known: ref = max(ref, known)
    rng = np.random.default_rng(5); out = dict(name=name, n=n, ref=float(ref), known=known, t_sa=t_sa, bsb={})
    for steps in (1000, 3000, 10000, 30000):
        cuts = []; t0 = time.perf_counter()
        for b in range(4):
            X = sb_batch(a, J, steps, 8, rng, False)
            cuts += [cut(Wm, 2.0 * x - 1.0) for x in X]
        out["bsb"][str(steps)] = dict(best=float(max(cuts)), median=float(np.median(cuts)), t=(time.perf_counter() - t0) / 4,
                                      gap_best=float((1 - max(cuts) / ref) * 100), gap_median=float((1 - np.median(cuts) / ref) * 100))
    with open("results/gset_check.jsonl", "a") as fh: fh.write(json.dumps(out) + "\n")
    print(name, "ref", ref, {k: (round(v["gap_best"], 2), round(v["gap_median"], 2)) for k, v in out["bsb"].items()}, flush=True)


if __name__ == "__main__":
    for k in sorted(GSET_BEST):
        if os.path.exists(f"gset/G{k}"): run(f"G{k}", read_gset(f"gset/G{k}"), GSET_BEST[k])
    if "--gset-only" not in sys.argv:
        for s in range(2): run(f"ER800-{s}", er(800, 0.06, s))
        for s in range(2): run(f"ER2000-{s}", er(2000, 0.01, s))
        for s in range(2): run(f"TOR800-{s}", torus(20, 40, s))

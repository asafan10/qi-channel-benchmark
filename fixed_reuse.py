"""The reuse baseline P10 takes, per instance, the best of the candidate lattice patterns (an oracle choice). Here one
pattern is fixed in advance for every instance: the candidate with the lowest cost on the unjittered lattice of the same size (first in
enumeration order on ties). -> results/fixed_reuse.json"""
import json, math, itertools, numpy as np
from exp_common import *
from refs import build_ref, sdp_bounds


def patterns(n, K):
    cols = int(math.ceil(math.sqrt(n * 1.15))); idx = np.arange(n); c = idx % cols; r = idx // cols; q = c - (r // 2); pats = []
    if K == 3:
        for a, b in itertools.product(range(3), range(3)):
            if (a, b) != (0, 0): pats.append(((a * q + b * r) % 3).astype(np.int64))
    else:
        for a, b, a2, b2 in itertools.product(range(2), repeat=4):
            p = 2 * ((a * q + b * r) % 2) + ((a2 * q + b2 * r) % 2)
            if len(set(p.tolist())) == 4: pats.append(p.astype(np.int64))
    return pats


def lattice_W(n):
    cols = int(math.ceil(math.sqrt(n * 1.15))); pts = []; r = 0
    while len(pts) < n:
        for c in range(cols): pts.append((c + 0.5 * (r % 2), r * math.sqrt(3) / 2))
        r += 1
    pts = np.array(pts[:n]); d = np.linalg.norm(pts[:, None] - pts[None], axis=2); np.fill_diagonal(d, np.inf)
    W = np.where(d <= 2.6, np.maximum(d, 0.5) ** -3.5, 0.0); return W


REF = build_ref(); LB = sdp_bounds(); out = {}
for K, n, seeds in ((3, 200, range(30)), (3, 500, range(30)), (4, 200, range(30))):
    P = patterns(n, K); W0 = lattice_W(n); c0 = csr_from_W(W0)
    fix = int(np.argmin([interference(p, *c0) for p in P]))
    best_g, fix_g, bb, fb = [], [], [], []
    for s in seeds:
        W = make_instance(n, s, K); csr = csr_from_W(W); costs = [interference(p, *csr) for p in P]
        best_g.append((min(costs) / REF[(K, n, s)] - 1) * 100); fix_g.append((costs[fix] / REF[(K, n, s)] - 1) * 100)
        if (K, n, s) in LB and LB[(K, n, s)]: bb.append((min(costs) / LB[(K, n, s)] - 1) * 100); fb.append((costs[fix] / LB[(K, n, s)] - 1) * 100)
    out[f"{K}_{n}"] = dict(fixed_index=fix, best_mean=float(np.mean(best_g)), fixed_mean=float(np.mean(fix_g)), diff_max=float(max(f - b for f, b in zip(fix_g, best_g))),
                           same=int(sum(abs(f - b) < 1e-9 for f, b in zip(fix_g, best_g))), n=len(best_g),
                           bound_best=float(np.mean(bb)) if bb else None, bound_fixed=float(np.mean(fb)) if fb else None, n_bound=len(bb))
json.dump(out, open("results/fixed_reuse.json", "w"), indent=1); print(json.dumps(out, indent=1))
from bench_jitter import make_layout
P = patterns(200, 3); fix = int(np.argmin([interference(p, *csr_from_W(lattice_W(200))) for p in P])); d = []
for s in range(30):
    W = make_layout(200, s, 3, 0.3); csr = csr_from_W(W); costs = [interference(p, *csr) for p in P]; d.append((costs[fix] / min(costs) - 1) * 100)
out["irregular_0.3"] = dict(same=int(sum(x < 1e-9 for x in d)), n=30, diff_max=float(max(d)), diff_median=float(np.median(d)))
json.dump(out, open("results/fixed_reuse.json", "w"), indent=1); print(out["irregular_0.3"])

"""Decoder-only floor with several scan orders, and weighted DSATUR.
For every Table XIV instance: sequential greedy (assign each cell the channel with least interference to the cells
already assigned) in (i) index order (= repair of the empty vector), (ii) reverse index order, (iii) decreasing weighted
degree (= the paper's Greedy), (iv) best of 100 uniformly random orders; plus weighted DSATUR (next cell = most distinct
channels among assigned neighbours, ties by weighted degree). -> results/decoder_floor.json"""
import json, numpy as np
from exp_common import setup, decode_eval, interference
from refs import build_ref


def order_greedy(order, K, ip, ix, dt_):
    n = len(order); ch = -np.ones(n, dtype=np.int64)
    for i in order:
        cost = np.zeros(K)
        for k in range(ip[i], ip[i + 1]):
            if ch[ix[k]] >= 0: cost[ch[ix[k]]] += dt_[k]
        ch[i] = int(np.argmin(cost))
    return ch


def dsatur(K, ip, ix, dt_):
    n = len(ip) - 1; ch = -np.ones(n, dtype=np.int64)
    deg = np.array([dt_[ip[i]:ip[i + 1]].sum() for i in range(n)])
    seen = [set() for _ in range(n)]
    for _ in range(n):
        best = None
        for i in range(n):
            if ch[i] >= 0: continue
            key = (len(seen[i]), deg[i])
            if best is None or key > best[0]: best = (key, i)
        i = best[1]; cost = np.zeros(K)
        for k in range(ip[i], ip[i + 1]):
            if ch[ix[k]] >= 0: cost[ch[ix[k]]] += dt_[k]
        ch[i] = int(np.argmin(cost))
        for k in range(ip[i], ip[i + 1]): seen[ix[k]].add(ch[i])
    return ch


if __name__ == "__main__":
    ref = build_ref(); out = {}
    for K, ns in ((3, (20, 50, 100, 200, 500)), (4, (50, 100, 200))):
        for n in ns:
            for s in range(5):
                W, csr, a, J, lam, jj = setup(n, s, K); ip, ix, dt_ = csr; r = ref[(K, n, s)]
                g = lambda ch: (interference(ch, *csr) / r - 1) * 100
                f0 = (decode_eval(np.zeros(n * K, dtype=np.int64), K, csr)[0] / r - 1) * 100
                idx = g(order_greedy(np.arange(n), K, *csr)); rev = g(order_greedy(np.arange(n)[::-1], K, *csr))
                deg = np.array([dt_[ip[i]:ip[i + 1]].sum() for i in range(n)])
                dg = g(order_greedy(np.argsort(-deg, kind="stable"), K, *csr))
                rng = np.random.default_rng(1000 * n + s); rnd = [g(order_greedy(rng.permutation(n), K, *csr)) for _ in range(100)]
                ds = g(dsatur(K, *csr))
                out[f"{K}_{n}_{s}"] = dict(repair_empty=f0, index=idx, reverse=rev, degree=dg, random_best=min(rnd),
                                           random_median=float(np.median(rnd)), dsatur=ds)
            v = [out[f"{K}_{n}_{s}"] for s in range(5)]
            print(K, n, {k: round(float(np.mean([x[k] for x in v])), 1) for k in v[0]}, flush=True)
    json.dump(out, open("results/decoder_floor.json", "w"), indent=1)

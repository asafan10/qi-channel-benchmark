"""Breakout local search (BLS) [Benlic and Hao 2013, published for Max-Cut; adapted here to K channels] as a state-of-the-art
classical baseline at matched wall-clock time. Moves change the channel of one cell; the local search is steepest descent over all
(cell, channel) moves; a local optimum is perturbed by L moves that are either directed (best non-tabu move, tabu tenure drawn from
[3, n/10], aspiration by the best cost) or random. L starts at L0 = max(2, 0.01 n) and grows by one each time the search returns to
the previous local optimum; the probability of a directed perturbation is max(exp(-w/T), P0) with w the number of consecutive
non-improving local optima, T = 1000 and P0 = 0.8, and after T non-improving optima a random perturbation of Lmax = 0.1 n moves is
applied. These are the published default settings; we did not tune them.
Runs are anytime: the best cost found by 1, 3, 10, 30 and 100 s of wall-clock time is recorded in one 100 s run from a random plan.
 T7   the instances of Table 7 (K = 3: 200, 500 cells; K = 4: 200 cells; seeds 0-4), three runs each  -> results/bls_T7.jsonl
 T12  the 30 irregular instances of Table 12 (jitter 0.3, K = 3, 200 cells), two runs each         -> results/bls_T12.jsonl
Usage: python3 bench_bls.py T7|T12 a-b [K n]"""
import sys, os, json, math, time, numpy as np
from numba import njit
from exp_common import *

BUD = [1, 3, 10, 30, 100]


@njit(cache=True)
def _move(i, c, ch, D, ip, ix, wd):
    o = ch[i]
    for k in range(ip[i], ip[i + 1]):
        j = ix[k]; D[j, o] -= wd[k]; D[j, c] += wd[k]
    ch[i] = c
    return D[i, c] - D[i, o]


@njit(cache=True)
def bls_chunk(ch, D, tabu, bch, st, ip, ix, wd, K, n_opt, seed):
    """Runs n_opt local-search/perturbation cycles. st = [f, fbest, L, w, it, fprev, L0, Lmax, T, P0]."""
    np.random.seed(seed)
    n = ch.shape[0]
    f, fbest, L, w, it, fprev, L0, Lmax, T, P0 = st[0], st[1], st[2], st[3], st[4], st[5], st[6], st[7], st[8], st[9]
    for _ in range(n_opt):
        while True:                                   # steepest descent
            bd = -1e-12; bi = -1; bc = -1
            for i in range(n):
                ci = ch[i]; di = D[i, ci]
                for c in range(K):
                    if c != ci and D[i, c] - di < bd: bd = D[i, c] - di; bi = i; bc = c
            if bi < 0: break
            f += _move(bi, bc, ch, D, ip, ix, wd)
        if f < fbest - 1e-12:
            fbest = f; w = 0.0
            for i in range(n): bch[i] = ch[i]
        else:
            w += 1.0
        if abs(f - fprev) < 1e-12: L += 1.0
        else: L = L0
        fprev = f
        if w > T:                                     # strong perturbation
            for _m in range(int(Lmax)):
                i = np.random.randint(n); c = np.random.randint(K - 1)
                if c >= ch[i]: c += 1
                f += _move(i, c, ch, D, ip, ix, wd)
            w = 0.0
        else:
            P = max(math.exp(-w / T), P0); directed = np.random.random() < P
            for _m in range(int(L)):
                it += 1.0
                if directed:
                    bd = 1e300; bi = -1; bc = -1
                    for i in range(n):
                        ci = ch[i]; di = D[i, ci]
                        for c in range(K):
                            if c == ci: continue
                            d = D[i, c] - di
                            if (tabu[i, c] <= it or f + d < fbest - 1e-12) and d < bd: bd = d; bi = i; bc = c
                    if bi < 0: continue
                    o = ch[bi]; f += _move(bi, bc, ch, D, ip, ix, wd)
                    tabu[bi, o] = it + 3 + np.random.randint(max(1, n // 10 - 2))
                else:
                    i = np.random.randint(n); c = np.random.randint(K - 1)
                    if c >= ch[i]: c += 1
                    f += _move(i, c, ch, D, ip, ix, wd)
    st[0], st[1], st[2], st[3], st[4], st[5] = f, fbest, L, w, it, fprev
    return fbest


def bls_run(W, csr, K, seed, tmax=100.0):
    ip, ix, wd = csr; n = W.shape[0]; rng = np.random.default_rng(seed)
    ch = rng.integers(0, K, n).astype(np.int64); D = np.zeros((n, K))
    for i in range(n):
        for k in range(ip[i], ip[i + 1]): D[i, ch[ix[k]]] += wd[k]
    f = 0.5 * sum(D[i, ch[i]] for i in range(n)); bch = ch.copy(); tabu = np.zeros((n, K))
    L0 = max(2.0, 0.01 * n); st = np.array([f, f, L0, 0.0, 0.0, -1.0, L0, max(L0, 0.1 * n), 1000.0, 0.8])
    t0 = time.perf_counter(); rec = {}; b = 0; k = 0
    while True:
        bls_chunk(ch, D, tabu, bch, st, ip, ix, wd, K, 20, int(rng.integers(1, 2**31 - 1))); k += 1
        el = time.perf_counter() - t0
        while b < len(BUD) and el >= BUD[b]:
            rec[str(BUD[b])] = float(interference(bch, ip, ix, wd)); b += 1
        if b == len(BUD): break
    assert abs(rec[str(BUD[-1])] - st[1]) < 1e-6 * max(1.0, st[1])
    return rec


def warm_bls():
    W = make_instance(20, 0, 3); csr = csr_from_W(W); ip, ix, wd = csr; n = 20
    ch = np.zeros(n, np.int64); D = np.zeros((n, 3))
    for i in range(n):
        for k in range(ip[i], ip[i + 1]): D[i, 0] += wd[k]
    bls_chunk(ch, D, np.zeros((n, 3)), ch.copy(), np.array([1.0, 1.0, 2, 0, 0, -1, 2, 2, 1000, 0.8]), ip, ix, wd, 3, 2, 1)


if __name__ == "__main__":
    part = sys.argv[1]; a_, b_ = map(int, sys.argv[2].split("-")); warm(3); warm_bls()
    if part == "T7":
        K, n = int(sys.argv[3]), int(sys.argv[4]); jl = f"results/bls_T7_K{K}_{n}.jsonl"
        done = {json.loads(l)["seed"] for l in open(jl)} if os.path.exists(jl) else set()
        for seed in range(a_, b_):
            if seed in done: continue
            W = make_instance(n, seed, K); csr = csr_from_W(W)
            runs = [bls_run(W, csr, K, 23000 + 100 * seed + r) for r in range(3)]
            append_jsonl(jl, dict(K=K, n=n, seed=seed, runs=runs)); print("T7", K, n, seed, [r["100"] for r in runs], flush=True)
    else:
        from bench_jitter import make_layout
        jl = "results/bls_T12.jsonl"; done = {json.loads(l)["seed"] for l in open(jl)} if os.path.exists(jl) else set()
        for seed in range(a_, b_):
            if seed in done: continue
            W = make_layout(200, seed, 3, 0.3); csr = csr_from_W(W)
            runs = [bls_run(W, csr, 3, 24000 + 100 * seed + r) for r in range(2)]
            append_jsonl(jl, dict(K=3, n=200, seed=seed, runs=runs)); print("T12", seed, [r["100"] for r in runs], flush=True)

"""Shared helpers for the additional experiments (A-D). Imports the frozen solvers from bench_ran_qubo."""
import json, math, time, os
import numpy as np
from numba import njit
from bench_ran_qubo import *  # make_instance, csr_from_W, build_qubo, sa_native, tabu, sb_batch, repair, interference, decode_eval, qubo_energy, CFG


@njit(cache=True)
def descent(ch0, K, ip, iy, dt_):
    """Steepest single-cell recolouring descent in native space until a local minimum (1-opt).
    Returns (assignment, cost, number of moves)."""
    n = len(ch0); ch = ch0.copy()
    cost = np.zeros((n, K))
    for i in range(n):
        for k in range(ip[i], ip[i + 1]):
            cost[i, ch[iy[k]]] += dt_[k]
    E = interference(ch, ip, iy, dt_); moves = 0
    while True:
        bd = -1e-12; bi = -1; bc = -1
        for i in range(n):
            ci = cost[i, ch[i]]
            for c in range(K):
                if c != ch[i]:
                    d = cost[i, c] - ci
                    if d < bd:
                        bd = d; bi = i; bc = c
        if bi < 0:
            break
        old = ch[bi]
        for k in range(ip[bi], ip[bi + 1]):
            j = iy[k]; cost[j, old] -= dt_[k]; cost[j, bc] += dt_[k]
        ch[bi] = bc; E += bd; moves += 1
    return ch, E, moves


def repair_all(Xb, K, csr):
    """Repair every replica of a bSB batch; returns list of (assignment, cost)."""
    out = []
    for x in Xb:
        ch, _ = repair(x.astype(np.int64), K, *csr)
        out.append((ch, interference(ch, *csr)))
    return out


def setup(n, seed, K, cfg=CFG):
    W = make_instance(n, seed, K); csr = csr_from_W(W)
    a, J, lam = build_qubo(W, K, cfg["lam_scale"])
    jp, jx, jd = J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64)
    return W, csr, a, J, lam, (jp, jx, jd)


def warm(K=3):
    W, csr, a, J, lam, jj = setup(12, 0, K)
    rng = np.random.default_rng(0)
    sa_native(rng.integers(0, K, 12), K, *csr, 2.0, 0.02, 5, 0)
    descent(rng.integers(0, K, 12), K, *csr)
    Xb = sb_batch(a, J, 20, 2, rng, False); qubo_energy(Xb[0], a, *jj); repair_all(Xb, K, csr)


def append_jsonl(path, rec):
    with open(path, "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def done_keys(path, keys=("n", "seed", "K")):
    s = set()
    if os.path.exists(path):
        for l in open(path):
            d = json.loads(l); s.add(tuple(d[k] for k in keys))
    return s

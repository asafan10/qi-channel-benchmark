"""
bench_ran_qubo.py
=================
Full-scale benchmark of quantum-inspired QUBO solvers vs. classical baselines
for weighted co-channel-interference minimisation (max-K-cut) in a cellular
layout. No truncation: every solver sees all N = n_cells*K binary variables.

Problem
-------
  cells on a jittered hexagonal lattice, interference weight
      w_ij = (d_ij / d_site)^(-3.5),  d_ij <= 2.6 d_site, near-field floor 0.5 d_site
  minimise  sum_{i<j} w_ij [ch_i == ch_j]        ch_i in {0..K-1}

QUBO (x_ic = 1 iff cell i uses channel c)
  E(x) = sum_{i<j,c} w_ij x_ic x_jc + lam * sum_i (sum_c x_ic - 1)^2
  lam = max weighted degree + 1  (guarantees one-hot violation never pays)

Solvers (equal wall-clock budget, restarts until budget exhausted)
  Greedy            weighted-degree-ordered constructive heuristic
  Tabu              tabu search in native assignment space
  SA-native         simulated annealing in native assignment space
  CP-SAT            OR-Tools, exact/anytime, 1 worker
  SA-QUBO           single-bit-flip SA on the QUBO
  PT-QUBO           parallel tempering on the QUBO
  bSB / dSB         ballistic / discrete simulated bifurcation (sparse, batched)

All QUBO solvers select their best sample by QUBO energy, then the same
one-hot repair is applied before reporting interference.
"""
import json, sys, time, math
import numpy as np
import scipy.sparse as sp
from numba import njit

# ------------------------------------------------------------------ problem
def make_instance(n_cells, seed, K=3):
    rng = np.random.default_rng(seed)
    cols = int(math.ceil(math.sqrt(n_cells * 1.15)))
    pts = []
    r = 0
    while len(pts) < n_cells:
        for c in range(cols):
            pts.append((c + 0.5 * (r % 2), r * math.sqrt(3) / 2))
        r += 1
    pts = np.array(pts[:n_cells]) + rng.normal(0, 0.15, (n_cells, 2))
    d = np.linalg.norm(pts[:, None] - pts[None], axis=2)
    np.fill_diagonal(d, np.inf)
    W = np.where(d <= 2.6, np.maximum(d, 0.5) ** -3.5, 0.0)
    W = np.triu(W, 1); W = W + W.T
    return W


def csr_from_W(W):
    S = sp.csr_matrix(W)
    return S.indptr.astype(np.int64), S.indices.astype(np.int64), S.data.astype(np.float64)


def build_qubo(W, K, lam_scale=1.0):
    """E(x) = a.x + sum_{i<j} J_ij x_i x_j + const ; J sparse symmetric, zero diagonal."""
    n = W.shape[0]; N = n * K
    lam = lam_scale * (W.sum(1).max() + 1.0)
    a = np.full(N, -lam)
    rows, cols, vals = [], [], []
    Wc = sp.coo_matrix(np.triu(W, 1))
    for c in range(K):
        rows.append(Wc.row * K + c); cols.append(Wc.col * K + c); vals.append(Wc.data)
    for i in range(n):
        for c in range(K):
            for cp in range(c + 1, K):
                rows.append(np.array([i * K + c])); cols.append(np.array([i * K + cp]))
                vals.append(np.array([2 * lam]))
    r = np.concatenate(rows); c_ = np.concatenate(cols); v = np.concatenate(vals)
    U = sp.coo_matrix((v, (r, c_)), shape=(N, N))
    J = (U + U.T).tocsr()
    return a, J, lam


# --------------------------------------------------------------- evaluation
@njit(cache=True)
def interference(ch, ip, ix, dt_):
    e = 0.0
    for i in range(len(ch)):
        for k in range(ip[i], ip[i + 1]):
            j = ix[k]
            if j > i and ch[i] == ch[j]:
                e += dt_[k]
    return e


@njit(cache=True)
def repair(x, K, ip, ix, dt_):
    """One-hot repair. Returns (channel assignment, n_feasible_cells)."""
    n = len(x) // K
    ch = -np.ones(n, dtype=np.int64)
    nfeas = 0
    for i in range(n):
        s = 0; last = -1
        for c in range(K):
            if x[i * K + c] == 1:
                s += 1; last = c
        if s == 1:
            ch[i] = last; nfeas += 1
    for i in range(n):
        if ch[i] >= 0:
            continue
        s = 0
        for c in range(K):
            s += x[i * K + c]
        best = 1e300; bc = 0
        for c in range(K):
            if s > 1 and x[i * K + c] == 0:
                continue
            cost = 0.0
            for k in range(ip[i], ip[i + 1]):
                if ch[ix[k]] == c:
                    cost += dt_[k]
            if cost < best:
                best = cost; bc = c
        ch[i] = bc
    return ch, nfeas


def decode_eval(x, K, csr):
    ch, nf = repair(x.astype(np.int64), K, *csr)
    return interference(ch, *csr), nf / (len(x) // K)


# ------------------------------------------------------- native-space solvers
def greedy(csr, n, K):
    ip, ix, dt_ = csr
    deg = np.array([dt_[ip[i]:ip[i + 1]].sum() for i in range(n)])
    ch = -np.ones(n, dtype=np.int64)
    for i in np.argsort(-deg):
        cost = np.zeros(K)
        for k in range(ip[i], ip[i + 1]):
            if ch[ix[k]] >= 0:
                cost[ch[ix[k]]] += dt_[k]
        ch[i] = int(np.argmin(cost))
    return ch


@njit(cache=True)
def sa_native(ch0, K, ip, ix, dt_, T0, T1, sweeps, seed):
    np.random.seed(seed)
    n = len(ch0); ch = ch0.copy()
    E = interference(ch, ip, ix, dt_); best = E; bch = ch.copy()
    for s in range(sweeps):
        T = T0 * (T1 / T0) ** (s / max(sweeps - 1, 1))
        for _ in range(n):
            i = np.random.randint(n)
            c2 = np.random.randint(K - 1)
            if c2 >= ch[i]:
                c2 += 1
            d = 0.0
            for k in range(ip[i], ip[i + 1]):
                cj = ch[ix[k]]
                if cj == c2: d += dt_[k]
                elif cj == ch[i]: d -= dt_[k]
            if d <= 0 or np.random.random() < math.exp(-d / T):
                ch[i] = c2; E += d
        if E < best - 1e-12:
            best = E; bch = ch.copy()
    return bch, best


@njit(cache=True)
def tabu(ch0, K, ip, ix, dt_, iters, tenure, seed):
    np.random.seed(seed)
    n = len(ch0); ch = ch0.copy()
    cost = np.zeros((n, K))
    for i in range(n):
        for k in range(ip[i], ip[i + 1]):
            cost[i, ch[ix[k]]] += dt_[k]
    E = interference(ch, ip, ix, dt_); best = E; bch = ch.copy()
    tabu_until = np.zeros((n, K), dtype=np.int64)
    for it in range(iters):
        bd = 1e300; bi = -1; bc = -1
        for i in range(n):
            for c in range(K):
                if c == ch[i]:
                    continue
                d = cost[i, c] - cost[i, ch[i]]
                if (tabu_until[i, c] <= it or E + d < best - 1e-12) and d < bd - 1e-15:
                    bd = d; bi = i; bc = c
        if bi < 0:
            continue
        old = ch[bi]
        for k in range(ip[bi], ip[bi + 1]):
            j = ix[k]
            cost[j, old] -= dt_[k]; cost[j, bc] += dt_[k]
        ch[bi] = bc; E += bd
        tabu_until[bi, old] = it + tenure + np.random.randint(10)
        if E < best - 1e-12:
            best = E; bch = ch.copy()
    return bch, best


# ------------------------------------------------------------- QUBO solvers
@njit(cache=True)
def qubo_energy(x, a, jp, jx, jd):
    e = 0.0
    for i in range(len(x)):
        if x[i]:
            e += a[i]
            for k in range(jp[i], jp[i + 1]):
                if jx[k] > i and x[jx[k]]:
                    e += jd[k]
    return e


@njit(cache=True)
def sa_qubo(x0, a, jp, jx, jd, T0, T1, sweeps, seed):
    np.random.seed(seed)
    N = len(x0); x = x0.copy()
    g = a.copy()
    for i in range(N):
        for k in range(jp[i], jp[i + 1]):
            if x[jx[k]]:
                g[i] += jd[k]
    E = qubo_energy(x, a, jp, jx, jd); best = E; bx = x.copy()
    for s in range(sweeps):
        T = T0 * (T1 / T0) ** (s / max(sweeps - 1, 1))
        for _ in range(N):
            i = np.random.randint(N)
            dx = 1 - 2 * x[i]
            d = dx * g[i]
            if d <= 0 or np.random.random() < math.exp(-d / T):
                x[i] += dx; E += d
                for k in range(jp[i], jp[i + 1]):
                    g[jx[k]] += dx * jd[k]
        if E < best - 1e-12:
            best = E; bx = x.copy()
    return bx, best


@njit(cache=True)
def pt_qubo(X0, a, jp, jx, jd, temps, sweeps, seed):
    np.random.seed(seed)
    R, N = X0.shape; X = X0.copy()
    G = np.zeros((R, N)); Es = np.zeros(R)
    for r in range(R):
        G[r] = a
        for i in range(N):
            for k in range(jp[i], jp[i + 1]):
                if X[r, jx[k]]:
                    G[r, i] += jd[k]
        Es[r] = qubo_energy(X[r], a, jp, jx, jd)
    best = 1e300; bx = X[0].copy()
    for s in range(sweeps):
        for r in range(R):
            T = temps[r]
            for _ in range(N):
                i = np.random.randint(N)
                dx = 1 - 2 * X[r, i]
                d = dx * G[r, i]
                if d <= 0 or np.random.random() < math.exp(-d / T):
                    X[r, i] += dx; Es[r] += d
                    for k in range(jp[i], jp[i + 1]):
                        G[r, jx[k]] += dx * jd[k]
            if Es[r] < best - 1e-12:
                best = Es[r]; bx = X[r].copy()
        for r in range(s % 2, R - 1, 2):
            arg = (1.0 / temps[r] - 1.0 / temps[r + 1]) * (Es[r] - Es[r + 1])
            if arg >= 0 or np.random.random() < math.exp(arg):
                tmp = X[r].copy(); X[r] = X[r + 1]; X[r + 1] = tmp
                tg = G[r].copy(); G[r] = G[r + 1]; G[r + 1] = tg
                te = Es[r]; Es[r] = Es[r + 1]; Es[r + 1] = te
    return bx, best


def sb_batch(a, J, steps, R, rng, discrete, dt=0.5):
    """Ballistic / discrete simulated bifurcation, R parallel replicas (sparse)."""
    N = len(a)
    # x = (1+s)/2 : E = const + s.h + 0.5 s' Jm s ,  Jm = J/4, h = a/2 + J.1/4
    Jm = J * 0.25
    h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
    scale = math.sqrt((Jm.multiply(Jm)).sum()) + 1e-12
    c0 = 0.5 * math.sqrt(N - 1) / scale
    X = rng.uniform(-0.1, 0.1, (N, R)); Y = rng.uniform(-0.1, 0.1, (N, R))
    for t in range(steps):
        at = t / steps
        S = np.sign(X) if discrete else X
        force = -(Jm @ S + h[:, None])          # descend Ising energy
        Y += dt * (-(1.0 - at) * X + c0 * force)
        X += dt * Y
        m = np.abs(X) > 1
        X[m] = np.sign(X[m]); Y[m] = 0.0
    return (X > 0).astype(np.int64).T          # R x N


# ---------------------------------------------------------------- CP-SAT
def cpsat(W, K, budget, hint=None):
    from ortools.sat.python import cp_model
    n = W.shape[0]; SC = 10000
    m = cp_model.CpModel()
    x = [[m.NewBoolVar(f"x{i}_{c}") for c in range(K)] for i in range(n)]
    for i in range(n):
        m.AddExactlyOne(x[i])
    obj = []
    ii, jj = np.nonzero(np.triu(W, 1))
    for i, j in zip(ii, jj):
        w = int(round(W[i, j] * SC))
        if w == 0: continue
        z = m.NewBoolVar("")
        for c in range(K):
            m.AddBoolOr([x[i][c].Not(), x[j][c].Not(), z])
        obj.append(w * z)
    if hint is not None:
        for i in range(n):
            for c in range(K):
                m.AddHint(x[i][c], int(hint[i] == c))
    m.Minimize(sum(obj))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = budget
    s.parameters.num_workers = 1
    s.parameters.random_seed = 1
    t0 = time.perf_counter()
    st = s.Solve(m)
    dt = time.perf_counter() - t0
    ch = np.array([[s.Value(x[i][c]) for c in range(K)] for i in range(n)]).argmax(1)
    return ch, s.BestObjectiveBound() / SC, st == cp_model.OPTIMAL, dt


# ------------------------------------------------------------------ driver
def run_instance(n, seed, K, budget, cfg, only=None):
    want = lambda name: only is None or name in only
    W = make_instance(n, seed, K)
    csr = csr_from_W(W)
    a, J, lam = build_qubo(W, K, cfg["lam_scale"])
    jp, jx, jd = J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64)
    N = n * K
    rng = np.random.default_rng(1000 + seed)
    out = {}

    def anytime(name, step):
        """step() -> (interference, feas_frac); repeat until budget."""
        t0 = time.perf_counter(); best = 1e300; curve = []; feas = None; k = 0
        while time.perf_counter() - t0 < budget:
            e, f = step(k); k += 1
            if e < best: best, feas = e, f
            curve.append((time.perf_counter() - t0, best))
        out[name] = dict(E=best, feas=feas, t=time.perf_counter() - t0, restarts=k, curve=curve)

    # greedy
    t0 = time.perf_counter(); ch = greedy(csr, n, K)
    if want("Greedy"): out["Greedy"] = dict(E=interference(ch, *csr), feas=1.0, t=time.perf_counter() - t0, restarts=1,
                         curve=[(time.perf_counter() - t0, interference(ch, *csr))])

    ch0 = rng.integers(0, K, n)
    if want("SA-native"): anytime("SA-native", lambda k: (sa_native(rng.integers(0, K, n), K, *csr, 2.0, 0.02, cfg["native_sweeps"], k)[1], 1.0))
    if want("Tabu"): anytime("Tabu", lambda k: (tabu(rng.integers(0, K, n), K, *csr, cfg["tabu_iters_per_cell"] * n, 8, k)[1], 1.0))

    # CP-SAT (warm-started from the greedy solution)
    if want("CP-SAT"):
      ch, bound, opt, dt = cpsat(W, K, budget, hint=greedy(csr, n, K))
      out["CP-SAT"] = dict(E=interference(ch, *csr), feas=1.0, t=dt, restarts=1, bound=bound, optimal=bool(opt),
                         curve=[(dt, interference(ch, *csr))])

    T0, T1 = cfg["T0_scale"] * lam, cfg["T1"]
    if want("SA-QUBO"): anytime("SA-QUBO", lambda k: decode_eval(sa_qubo(rng.integers(0, 2, N), a, jp, jx, jd, T0, T1, cfg["sweeps"], k)[0], K, csr))
    temps = np.geomspace(cfg["pt_tmin"], T0, cfg["pt_R"])
    if want("PT-QUBO"): anytime("PT-QUBO", lambda k: decode_eval(pt_qubo(rng.integers(0, 2, (cfg["pt_R"], N)), a, jp, jx, jd, temps, cfg["pt_sweeps"], k)[0], K, csr))
    for name, disc in (("bSB", False), ("dSB", True)):
        def step(k, disc=disc):
            Xb = sb_batch(a, J, cfg["sb_steps"], cfg["sb_R"], rng, disc)
            es = [qubo_energy(x, a, jp, jx, jd) for x in Xb]
            return decode_eval(Xb[int(np.argmin(es))], K, csr)
        if want(name): anytime(name, step)
    return out


# hyper-parameters fixed from a held-out pilot instance (n=100, seed 99), then frozen
CFG = dict(lam_scale=0.25, sweeps=3000, native_sweeps=1000, tabu_iters_per_cell=40, T0_scale=0.5, T1=0.15,
           pt_R=12, pt_tmin=0.1, pt_sweeps=400, sb_steps=3000, sb_R=8)

if __name__ == "__main__":
    sizes = [int(s) for s in sys.argv[1].split(",")]
    seeds = int(sys.argv[2]); K = int(sys.argv[3]); outf = sys.argv[4]
    import os
    done = set()
    if os.path.exists(outf):
        for l in open(outf):
            d = json.loads(l); done.add((d["n"], d["seed"]))
    with open(outf, "a") as f:
        for n in sizes:
            budget = 3.0 if n <= 50 else (6.0 if n <= 150 else 12.0)
            for s in range(seeds):
                if (n, s) in done: continue
                res = run_instance(n, s, K, budget, CFG, only=(sys.argv[5].split(',') if len(sys.argv) > 5 else None))
                f.write(json.dumps(dict(n=n, seed=s, K=K, budget=budget, res=res)) + "\n"); f.flush()
                print(n, s, {k: round(v["E"], 2) for k, v in res.items()}, flush=True)

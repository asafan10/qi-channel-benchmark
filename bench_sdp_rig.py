"""Rigorous SDP bounds (SCS is a first-order solver, so its objective is not a certificate).
Re-solve the relaxation of bench_sdp.py, keep the dual multipliers y (diag X = 1) and mu >= 0 (X_ij >= -1/(K-1) on edges),
and evaluate weak duality with an eigenvalue correction, which is valid for ANY y and mu >= 0:
   max <C,X> + const  <=  const + sum(y) + sum(mu)/(K-1) + n * max(0, -lambda_min(Z)),  Z = Diag(y) - sum mu_e E_e - C,
because <Z,X> >= lambda_min(Z) tr(X) = n lambda_min(Z) for every feasible X.  -> results/sdp_rig_K*.jsonl"""
from ortools.sat.python import cp_model  # noqa: F401
import os, sys, time, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, cvxpy as cp
from bench_ran_qubo import make_instance


def rig(W, K, eps=1e-6):
    n = W.shape[0]; Wu = np.triu(W, 1); Wtot = Wu.sum(); ii, jj = np.nonzero(Wu)
    X = cp.Variable((n, n), PSD=True)
    cons = [cp.diag(X) == 1, X[ii, jj] >= -1.0 / (K - 1)]
    prob = cp.Problem(cp.Maximize((K - 1) / K * cp.sum(cp.multiply(Wu, 1 - X))), cons)
    t0 = time.perf_counter(); prob.solve(solver="SCS", eps=eps, max_iters=200000); t = time.perf_counter() - t0
    const = (K - 1) / K * Wtot
    C = np.zeros((n, n)); C[ii, jj] = -(K - 1) / (2 * K) * Wu[ii, jj]; C = C + C.T
    y0 = np.asarray(cons[0].dual_value).ravel(); mu0 = np.asarray(cons[1].dual_value).ravel()
    best = None
    for sy in (1, -1):
        for sm in (1, -1):
            y = sy * y0; mu = np.clip(sm * mu0, 0, None)
            Z = np.diag(y) - C; M = np.zeros((n, n)); M[ii, jj] = mu / 2; M = M + M.T; Z = Z - M
            lmin = float(np.linalg.eigvalsh(Z)[0])
            ub = const + y.sum() + mu.sum() / (K - 1) + n * max(0.0, -lmin)
            if best is None or ub < best[0]: best = (ub, lmin, sy, sm)
    ub, lmin, sy, sm = best
    return dict(W=float(Wtot), sdp=float(prob.value), ub_cut_rig=float(ub), lb_f_rig=float(Wtot - ub), lmin=lmin,
                slack=float(ub - prob.value), status=prob.status, t=t)


if __name__ == "__main__":
    plan = [(3, n) for n in (20, 50, 100, 200)] + [(4, n) for n in (50, 100, 200)]
    out = "results/sdp_rig.jsonl"; done = set()
    if os.path.exists(out):
        for l in open(out):
            d = json.loads(l); done.add((d["K"], d["n"], d["seed"]))
    for K, n in plan:
        for s in range(5):
            if (K, n, s) in done: continue
            r = rig(make_instance(n, s, K), K)
            with open(out, "a") as fh: fh.write(json.dumps(dict(K=K, n=n, seed=s, **r)) + "\n")
            print(K, n, s, "lb_rig %.4f sdp %.4f slack %.2e lmin %.2e %.0fs" % (r["lb_f_rig"], r["W"] - r["sdp"], r["slack"], r["lmin"], r["t"]), flush=True)

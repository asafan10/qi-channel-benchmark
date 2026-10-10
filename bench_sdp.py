"""
bench_sdp.py -- Experiment B: certified lower bounds on the optimal interference f* from the max-K-cut SDP relaxation
(Frieze & Jerrum 1997):  max  (K-1)/K * sum_{i<j} w_ij (1 - X_ij)  s.t.  X PSD, X_ii = 1, X_ij >= -1/(K-1) for (i,j) in E.
(Imposing the inequality only on edges relaxes the original SDP further, so the bound stays valid.)
Since f(c) = W - cut(c), every feasible plan satisfies f >= W - SDP*. The bound is made conservative by using the larger of
the solver's primal and dual objectives. CP-SAT's bound was 0.00 after 60 s at 50 and 100 cells (bound_probe.py), so it is not used.
"""
from ortools.sat.python import cp_model  # noqa: F401  (import before cvxpy: avoids a HiGHS symbol clash)
import os, sys, time, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, cvxpy as cp
from bench_ran_qubo import make_instance


def sdp_bound(W, K, eps=1e-6):
    n = W.shape[0]; Wu = np.triu(W, 1); Wtot = Wu.sum(); ii, jj = np.nonzero(Wu)
    X = cp.Variable((n, n), PSD=True)
    cons = [cp.diag(X) == 1, X[ii, jj] >= -1.0 / (K - 1)]
    prob = cp.Problem(cp.Maximize((K - 1) / K * cp.sum(cp.multiply(Wu, 1 - X))), cons)
    t0 = time.perf_counter(); prob.solve(solver="SCS", eps=eps, max_iters=200000)
    info = prob.solver_stats.extra_stats.get("info", {}) if prob.solver_stats.extra_stats else {}
    # SCS minimises the negated objective: pobj, dobj refer to the minimisation; upper bound on max-cut = max(-pobj, -dobj, value)
    cand = [prob.value] + [-info[k] for k in ("pobj", "dobj") if k in info]
    ub_cut = max(cand)
    return dict(W=float(Wtot), sdp=float(prob.value), ub_cut=float(ub_cut), lb_f=float(Wtot - ub_cut), status=prob.status,
                t=time.perf_counter() - t0, iters=int(info.get("iter", -1)))


if __name__ == "__main__":
    K = int(sys.argv[1]); sizes = [int(x) for x in sys.argv[2].split(",")]; outf = sys.argv[3]
    done = set()
    if os.path.exists(outf):
        for l in open(outf):
            d = json.loads(l); done.add((d["n"], d["seed"]))
    for n in sizes:
        for seed in range(5):
            if (n, seed) in done: continue
            r = sdp_bound(make_instance(n, seed, K), K)
            with open(outf, "a") as fh: fh.write(json.dumps(dict(n=n, seed=seed, K=K, **r)) + "\n")
            print(K, n, seed, "lb_f %.3f (sdp %.3f, %s, %.0fs, %d it)" % (r["lb_f"], r["sdp"], r["status"], r["t"], r["iters"]), flush=True)

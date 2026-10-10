"""SDP rounding baseline: Frieze-Jerrum rounding of the max-k-cut relaxation of Section 3.4 (P9),
followed by native-space descent and by a short native SA polish (T 0.3 -> 0.02, as in P6).
K = 3, n in {20, 50, 100} (seeds 0-4).  Each rounding draws K Gaussian vectors g_c and assigns cell i to argmax_c <v_i, g_c>,
where V V^T = X is the SDP solution (SCS, eps 1e-5).  -> results/sdp_round.jsonl"""
from ortools.sat.python import cp_model  # noqa: F401
import os, sys, time, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, cvxpy as cp
from bench_ran_qubo import make_instance, csr_from_W, interference, sa_native
from exp_common import descent
import refs
SIGMA = float(os.environ.get("SIGMA", "0"))
if SIGMA:
    from bench_jitter import make_layout
    make_instance = lambda n, seed, K: make_layout(n, seed, K, SIGMA)

R_ROUND, N_POLISH, POL_SWEEPS = 1000, 20, 4000


def solve_sdp(W, K, eps=1e-5):
    n = W.shape[0]; Wu = np.triu(W, 1); ii, jj = np.nonzero(Wu)
    X = cp.Variable((n, n), PSD=True)
    prob = cp.Problem(cp.Maximize((K - 1) / K * cp.sum(cp.multiply(Wu, 1 - X))), [cp.diag(X) == 1, X[ii, jj] >= -1.0 / (K - 1)])
    t0 = time.perf_counter(); prob.solve(solver="SCS", eps=eps, max_iters=200000); return X.value, time.perf_counter() - t0


def main(K=3, sizes=(20, 50, 100), seeds=range(5), eps=1e-5):
    ref = refs.build_ref() if not SIGMA else {}; out = "results/sdp_round.jsonl" if not SIGMA else f"results/sdp_round_sigma{SIGMA}.jsonl"
    done = set()
    if os.path.exists(out):
        for l in open(out): d = json.loads(l); done.add((d["n"], d["seed"]))
    for n in sizes:
        for s in seeds:
            if (n, s) in done: continue
            W = make_instance(n, s, K); csr = csr_from_W(W); Xv, t_sdp = solve_sdp(W, K, eps)
            w, U = np.linalg.eigh((Xv + Xv.T) / 2); V = U * np.sqrt(np.clip(w, 0, None))
            rng = np.random.default_rng(1000 + s); t0 = time.perf_counter(); f_r, f_d, plans = [], [], []
            for _ in range(R_ROUND):
                ch = np.argmax(V @ rng.normal(size=(n, K)), axis=1).astype(np.int64)
                f_r.append(interference(ch, *csr)); chd, fd, _ = descent(ch, K, *csr); f_d.append(fd); plans.append(chd)
            t_round = time.perf_counter() - t0
            order = np.argsort(f_d)[:N_POLISH]; t0 = time.perf_counter()
            f_p = [sa_native(plans[k].copy(), K, *csr, 0.3, 0.02, POL_SWEEPS, 77 + int(k))[1] for k in order]
            t_pol = time.perf_counter() - t0
            fr = ref.get((K, n, s), 1.0); g = lambda v: float((v / fr - 1) * 100)
            rec = dict(K=K, n=n, seed=s, ref=fr, t_sdp=t_sdp, t_round=t_round, t_polish=t_pol,
                       round_best=g(min(f_r)), round_median=g(np.median(f_r)), desc_best=g(min(f_d)), desc_median=g(np.median(f_d)),
                       polish_best=g(min(f_p)), polish_median=g(np.median(f_p)), f_best=float(min(min(f_p), min(f_d))), f_round=float(min(f_r)), f_desc=float(min(f_d)), f_desc_med=float(np.median(f_d)), f_pol=float(min(f_p)), sigma=SIGMA, eps=eps, plan=[int(c) for c in plans[int(np.argmin(f_d))]] if min(f_d) <= min(f_p) else None)
            with open(out, "a") as fh: fh.write(json.dumps(rec) + "\n")
            print(n, s, {k: round(v, 3) for k, v in rec.items() if isinstance(v, float)}, flush=True)


if __name__ == "__main__":
    a = sys.argv[1:]
    main(sizes=tuple(int(x) for x in a[0].split(",")) if a else (20, 50, 100),
         seeds=[int(x) for x in a[1].split(",")] if len(a) > 1 else range(5), eps=float(a[2]) if len(a) > 2 else 1e-5)

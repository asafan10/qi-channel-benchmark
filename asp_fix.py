import json, time, random, collections, numpy as np
from bench_ran_qubo import *
import quantum_ran_complete as q0, qrc_fix1 as q1, qrc_fix2 as q2
ref = json.load(open("final_reference.json"))["K3"]; K = 3; rows = []
for n in (20, 50):
    for seed in range(5):
        W = make_instance(n, seed, K); csr = csr_from_W(W); N = n * K
        a, J, lam = build_qubo(W, K, CFG["lam_scale"]); Q = J.toarray() / 2.0; np.fill_diagonal(Q, a)
        for name, sv in (("DA as-given", q0.DA(seed=seed)), ("DA + dE fix", q1.DA(seed=seed)), ("DA + dE fix + swap fix", q2.DA(seed=seed)),
                         ("PA as-given", q0.PA(seed=seed)), ("PA + dE fix", q1.PA(seed=seed))):
            t0 = time.perf_counter(); x, e = sv.solve(Q, max_n=N); dt = time.perf_counter() - t0
            E, f = decode_eval(np.asarray(x).astype(np.int64), K, csr)
            true_e = float(np.asarray(x) @ Q @ np.asarray(x))
            rows.append(dict(n=n, seed=seed, solver=name, gap=100 * (E - ref[f"{n}_{seed}"]) / ref[f"{n}_{seed}"], feas=f, t=dt, reported_minus_true=float(e) - true_e))
    print("done", n, flush=True)
json.dump(rows, open("asp_fix.json", "w"))
g = collections.defaultdict(list)
for r in rows: g[(r["n"], r["solver"])].append(r)
print("| n | solver | gap % | raw one-hot feasibility | reported energy - true energy (mean) | time (s) |\n|---|---|---|---|---|---|")
for n in (20, 50):
    for s in ("DA as-given", "DA + dE fix", "DA + dE fix + swap fix", "PA as-given", "PA + dE fix"):
        v = g[(n, s)]; print(f"| {n} | {s} | {np.mean([r['gap'] for r in v]):.1f} | {np.mean([r['feas'] for r in v]):.2f} | {np.mean([r['reported_minus_true'] for r in v]):.1f} | {np.mean([r['t'] for r in v]):.1f} |")

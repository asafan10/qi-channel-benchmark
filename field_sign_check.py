"""Is the sign of the field in our bSB correct?  On small instances with a nonzero field (the one-hot QUBO,
rho = 0.25) the optimum is found by enumerating all K^n channel plans (Prop. 1 makes the Ising ground state one-hot).  Our bSB
(bench_ran_qubo.sb_batch, unchanged) is run, and a copy with the sign of the field flipped; both are decoded and repaired as in the paper.
Also checks that the Ising energy of Eq. (10) at the decoded spins equals the QUBO energy, and that it decreases along the run.
-> results/field_sign_check.json"""
import json, itertools, math, numpy as np
from exp_common import *
from bench_ran_qubo import sb_batch, qubo_energy


def sb_flipped(a, J, steps, R, rng):
    N = len(a); Jm = J * 0.25; h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
    c0 = 0.5 * math.sqrt(N - 1) / math.sqrt((Jm.multiply(Jm)).sum()); X = rng.uniform(-0.1, 0.1, (N, R)); Y = rng.uniform(-0.1, 0.1, (N, R))
    for t in range(steps):
        Y += 0.5 * (-(1.0 - t / steps) * X + c0 * (-(Jm @ X) + h[:, None])); X += 0.5 * Y      # field sign flipped
        m = np.abs(X) > 1; X[m] = np.sign(X[m]); Y[m] = 0.0
    return (X > 0).astype(np.int64).T


out = []
for n in (8, 10, 12):
    for seed in range(5):
        K = 3; W, csr, a, J, lam, jj = setup(n, seed, K)
        best = min(interference(np.array(c, dtype=np.int64), *csr) for c in itertools.product(range(K), repeat=n))
        res = {}
        for name, fn in (("ours", lambda r: sb_batch(a, J, 3000, 8, r, False)), ("field_flipped", lambda r: sb_flipped(a, J, 3000, 8, r))):
            hits, feas, gaps = 0, [], []
            for b in range(10):
                X = fn(np.random.default_rng(100 * seed + b))
                es = [qubo_energy(x, a, *jj) for x in X]; x = X[int(np.argmin(es))]
                ch, nf = repair(x.astype(np.int64), K, *csr); f = interference(ch, *csr)
                feas.append(float(np.mean(x.reshape(n, K).sum(1) == 1))); gaps.append((f / best - 1) * 100); hits += int(f <= best + 1e-9)
            res[name] = dict(optimum_hits_of_10=hits, raw_feasibility=float(np.mean(feas)), median_gap=float(np.median(gaps)))
        # Ising energy of Eq. (10) equals the QUBO energy on decoded states
        Jm = (J * 0.25).toarray(); h = a / 2 + np.asarray(J.sum(1)).ravel() / 4; x = sb_batch(a, J, 3000, 1, np.random.default_rng(seed), False)[0]
        s = 2.0 * x - 1.0; E_qubo = qubo_energy(x, a, *jj); E_ising = float(s @ h + 0.5 * s @ Jm @ s)
        res["energy_offset_check"] = float(E_qubo - E_ising)          # constant E0 (independent of the state)
        out.append(dict(n=n, seed=seed, optimum=float(best), **res)); print(n, seed, res, flush=True)
json.dump(out, open("results/field_sign_check.json", "w"), indent=1)

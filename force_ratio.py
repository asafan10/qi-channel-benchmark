"""Size of the penalty terms relative to the interference terms inside our bSB integrator (rho = 0.25).
Ising form of Section 3.5: E = s.h + 0.5 s' Jm s, Jm = J/4, h = a/2 + J.1/4. J splits into intra-cell penalty couplings (2 lambda) and
inter-cell interference couplings (w_ij); h splits into the penalty part lambda (K-2)/2 and the interference part delta_i/4.
Reports Frobenius-norm ratios and, along the pump of 3,000-step runs (steps 0, 750, 1,500, 2,250 and the final state), the median over spins
of |penalty force| / |interference force|, where the force on spin p is -c0 (Jm x + h)_p split into its penalty and interference parts.
-> results/force_ratio.json"""
import json, numpy as np, scipy.sparse as sp
import exp_common
from exp_common import *
from bench_jitter import make_layout
out = {}
for lay in ("regular", "jitter0.3"):
    if lay == "jitter0.3": exp_common.make_instance = lambda n, seed, K: make_layout(n, seed, K, 0.3)
    rows = []
    for seed in range(5):
        W, csr, a, J, lam, jj = exp_common.setup(200, seed, 3); K = 3; n = 200; N = n * K
        J = J.tocoo(); cell = lambda p: p // K
        pen = cell(J.row) == cell(J.col)
        Jp = sp.coo_matrix((J.data * pen, (J.row, J.col)), shape=(N, N)).tocsr() * 0.25
        Ji = sp.coo_matrix((J.data * ~pen, (J.row, J.col)), shape=(N, N)).tocsr() * 0.25
        hp = np.full(N, lam * (K - 2) / 2); hi = a / 2 + np.asarray(J.tocsr().sum(1)).ravel() / 4 - hp
        Jm = J.tocsr() * 0.25; h = hp + hi; c0 = 0.5 * np.sqrt(N - 1) / np.sqrt((Jm.multiply(Jm)).sum()); dt = 0.5; steps = 3000
        rng = np.random.default_rng(seed); Xs = rng.uniform(-0.1, 0.1, (N, 8)); Ys = rng.uniform(-0.1, 0.1, (N, 8)); snap = {}
        for t in range(steps + 1):
            if t in (0, 750, 1500, 2250, steps):
                fp = np.abs(Jp @ Xs + hp[:, None]); fi = np.abs(Ji @ Xs + hi[:, None]); snap[str(t)] = float(np.median(fp / np.maximum(fi, 1e-12)))
            if t == steps: break
            Ys += dt * (-(1.0 - t / steps) * Xs - c0 * (Jm @ Xs + h[:, None])); Xs += dt * Ys
            m_ = np.abs(Xs) > 1; Xs[m_] = np.sign(Xs[m_]); Ys[m_] = 0.0
        r_force = [snap[str(steps)]]
        rows.append(dict(seed=seed, lam=float(lam), frob_J=float(np.sqrt((Jp.multiply(Jp)).sum() / (Ji.multiply(Ji)).sum())),
                         field=float(np.linalg.norm(hp) / np.linalg.norm(hi)), force_median=float(np.median(r_force)), along_pump=snap))
    out[lay] = dict(rows=rows, **{k: float(np.median([r[k] for r in rows])) for k in ("lam", "frob_J", "field", "force_median")},
                    along_pump={t: float(np.median([r["along_pump"][t] for r in rows])) for t in rows[0]["along_pump"]})
    print(lay, {k: round(v, 3) for k, v in out[lay].items() if k != "rows"})
json.dump(out, open("results/force_ratio.json", "w"), indent=1)

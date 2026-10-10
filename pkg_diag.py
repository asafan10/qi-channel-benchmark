"""Diagnostic of the SB package collapse.
(1) Energy check: the package's own model, evaluated on our decoded states, against our Ising energy.
(2) Scale: the package normalises its coupling by the Frobenius norm of the matrix WITH the field embedded as an extra
    spin; we normalise by the couplings only.  Ratio of the two, and of field to coupling force on the auxiliary spin.
(3) Swap test at 200 cells: the package with its coupling scale replaced by ours, and our bSB with the package's scale.
-> results/pkg_diag.json"""
import json, math, warnings, numpy as np, torch, simulated_bifurcation as sb
from simulated_bifurcation.optimizer import simulated_bifurcation_optimizer as SBO
warnings.filterwarnings("ignore"); torch.set_num_threads(1)
from exp_common import *
from bench_ran_qubo import sb_batch

K = 3; out = {}


def pkg_inputs(a, J):
    Jm = (J * 0.25).toarray(); h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
    return Jm, h, torch.tensor(Jm / 2, dtype=torch.float32), torch.tensor(h, dtype=torch.float32)


def scales(Q, l):
    """Package tensor T (field embedded as an extra spin).  Package rule: 0.5 sqrt(N) / ||T||_F.
    Our rule applied to the same tensor: 0.5 sqrt(N-1) / ||T_couplings||_F.  Returns (ours, pkg, aux ratio, field/coupling Frobenius)."""
    T = sb.build_model(Q, l).to_ising(domain="spin").as_simulated_bifurcation_tensor().double().numpy(); N = T.shape[0] - 1
    C = T[:N, :N]; hf = T[:N, N]
    ours = 0.5 * math.sqrt(N - 1) / math.sqrt((C ** 2).sum()); pkg = 0.5 * math.sqrt(N) / math.sqrt((T ** 2).sum())
    return ours, pkg, float(np.abs(hf).sum() / np.abs(C).sum(1).max()), float(2 * (hf ** 2).sum() / (C ** 2).sum())


def sb_c0(a, J, steps, R, rng, mult):
    """Our bSB (bench_ran_qubo.sb_batch) with its coupling scale c0 multiplied by mult."""
    N = len(a); Jm = J * 0.25; h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
    c0 = mult * 0.5 * math.sqrt(N - 1) / math.sqrt((Jm.multiply(Jm)).sum()); dt = 0.5
    X = rng.uniform(-0.1, 0.1, (N, R)); Y = rng.uniform(-0.1, 0.1, (N, R))
    for t in range(steps):
        at = t / steps; Y += dt * (-(1.0 - at) * X - c0 * (Jm @ X + h[:, None])); X += dt * Y
        m = np.abs(X) > 1; X[m] = np.sign(X[m]); Y[m] = 0.0
    return (X > 0).astype(np.int64).T


_orig = SBO.SimulatedBifurcationOptimizer._SimulatedBifurcationOptimizer__init_quadratic_scale_parameter
def run_pkg(Q, l, csr, override=None, agents=128, ts=0.1, ps=0.003, steps=3333, seed=0):
    if override is not None:
        def f(self, matrix): self.quadratic_scale_parameter = torch.tensor(override, dtype=matrix.dtype)
        SBO.SimulatedBifurcationOptimizer._SimulatedBifurcationOptimizer__init_quadratic_scale_parameter = f
    try:
        sb.set_env(time_step=ts, pressure_slope=ps); torch.manual_seed(seed)
        S, E = sb.minimize(Q, l, domain="spin", agents=agents, max_steps=steps, mode="ballistic", heated=False,
                           verbose=False, best_only=False, early_stopping=False)
    finally:
        SBO.SimulatedBifurcationOptimizer._SimulatedBifurcationOptimizer__init_quadratic_scale_parameter = _orig; sb.reset_env()
    X = ((S.cpu().numpy() + 1) // 2).astype(np.int64)
    r = [decode_eval(x, K, csr) for x in X]; return float(min(f for f, _ in r)), float(np.mean([p for _, p in r]))


# (1) energy check
chk = []
for n in (50, 200):
    for s in range(3):
        W, csr, a, J, lam, jj = setup(n, s, K); Jm, h, Q, l = pkg_inputs(a, J)
        model = sb.build_model(Q, l)                       # x^T Q x + l.x  on the spin domain
        X = sb_batch(a, J, 3000, 4, np.random.default_rng(s), False)
        for x in X:
            sgn = 2.0 * x - 1.0
            ours = float(sgn @ h + 0.5 * sgn @ Jm @ sgn)
            theirs = float(model(torch.tensor(sgn, dtype=torch.float32)))
            chk.append(abs(ours - theirs) / max(1.0, abs(ours)))
out["energy_check_max_rel_diff"] = max(chk)
print("energy check, max relative difference:", max(chk))

# (2) scales
sc = {}
for n in (50, 100, 200, 500):
    W, csr, a, J, lam, jj = setup(n, 0, K); Jm, h, Q, l = pkg_inputs(a, J)
    o, p, aux_ratio, frob_ratio = scales(Q, l)
    sc[n] = dict(c0_ourrule=o, c0_pkg=p, pkg_over_ours=p / o, aux_force_over_max_spin_force=aux_ratio, field_frob_over_coupling_frob=frob_ratio)
    print(n, {k: round(v, 4) for k, v in sc[n].items()})
out["scales"] = sc

# (3) swap test at 200 cells, five instances: best decoded cost, mean raw feasibility
sw = {k: [] for k in ("pkg_default", "pkg_with_our_scale", "ours_default", "ours_with_pkg_scale")}
for s in range(5):
    W, csr, a, J, lam, jj = setup(200, s, K); Jm, h, Q, l = pkg_inputs(a, J)
    o, p, _, _ = scales(Q, l)
    sw["pkg_default"].append(run_pkg(Q, l, csr, seed=s))
    sw["pkg_with_our_scale"].append(run_pkg(Q, l, csr, override=o, seed=s))
    for key, mult in (("ours_default", 1.0), ("ours_with_pkg_scale", p / o)):
        Xb = sb_c0(a, J, 3000, 8, np.random.default_rng(s), mult)
        r = [decode_eval(x, K, csr) for x in Xb]; sw[key].append((float(min(f for f, _ in r)), float(np.mean([q for _, q in r]))))
    print(s, {k: (round(v[-1][0], 2), round(v[-1][1], 2)) for k, v in sw.items()}, flush=True)
out["swap200"] = sw
json.dump(out, open("results/pkg_diag.json", "w"), indent=1)

"""Minimal reproduction of the collapse of the simulated-bifurcation package (version 2.0.0) on the one-hot QUBO of this study
(Section 5.7 and Section S2 of Additional file 1). One 200-cell instance, K = 3, rho = 0.25; ballistic mode, eight agents, 3,000 steps,
time step 0.5 (our settings) and the package's default time step. Prints the fraction of cells with exactly one active channel (raw
feasibility) of the returned spins; our own bSB on the same Ising problem returns 1.00.
Usage: python3 pkg_collapse_repro.py"""
import numpy as np, torch, simulated_bifurcation as sb
from bench_ran_qubo import make_instance, build_qubo

K, n = 3, 200
W = make_instance(n, 0, K); a, J, lam = build_qubo(W, K, 0.25)               # E(x) = a.x + sum_{p<q} J_pq x_p x_q + const
Jm = (J * 0.25).toarray(); h = a / 2 + np.asarray(J.sum(1)).ravel() / 4       # Ising: E(s) = s.h + sum_{p<q} Jm_pq s_p s_q + const
Q = torch.tensor(Jm / 2, dtype=torch.float32); l = torch.tensor(h, dtype=torch.float32)   # package convention: s'Q s + l's
torch.manual_seed(0); torch.set_num_threads(1)
for label, env in (("time step 0.5", dict(time_step=0.5, pressure_slope=1.0 / (0.5 * 3000))), ("package defaults", {})):
    if env: sb.set_env(**env)
    s, e = sb.minimize(Q, l, domain="spin", agents=8, max_steps=3000, mode="ballistic", heated=False, verbose=False, best_only=True,
                       early_stopping=False)
    x = ((s.numpy().astype(int) + 1) // 2).reshape(n, K)
    print(f"{label}: raw feasibility {np.mean(x.sum(1) == 1):.2f}, active bits {int(x.sum())} of {n * K}")
    sb.reset_env()

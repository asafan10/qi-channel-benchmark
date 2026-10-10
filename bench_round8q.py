"""Planted 3-colourings at 50 and 100 cells, where the open-source SB package stays feasible (Section 4.8): an independent
implementation of bSB with fields on instances with known optimum 0. Our bSB (default), the package's ballistic mode with its
held-out settings (time step 0.1, pressure slope 0.003, 3,333 steps, 128 agents), and native SA. -> results/planted_small.jsonl"""
import time, numpy as np, torch, simulated_bifurcation as sb
from exp_common import *
from bench_round8 import sb_c0
torch.set_num_threads(1)
K = 3; out = "results/planted_small.jsonl"; warm(K)
for n in (50, 100, 200):
    for seed in range(5):
        W0 = make_instance(n, seed, K); col = np.random.default_rng(900 + seed).integers(0, K, n)
        W = W0 * (col[:, None] != col[None, :]); csr = csr_from_W(W)
        a, J, lam = build_qubo(W, K, 0.25); rng = np.random.default_rng(seed); res = {}
        res["bSB|3000"] = [float(decode_eval(x, K, csr)[0]) for b in range(4) for x in sb_c0(a, J, 3000, 8, rng, "default")]
        Jm = (J * 0.25).toarray(); h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
        Q = torch.tensor(Jm / 2, dtype=torch.float32); l = torch.tensor(h, dtype=torch.float32)
        sb.set_env(time_step=0.1, pressure_slope=0.003); pk = []; ph = []
        for k in range(3):
            torch.manual_seed(1000 * seed + k)
            S, E = sb.minimize(Q, l, domain="spin", agents=128, max_steps=3333, mode="ballistic", heated=False, verbose=False, best_only=True, early_stopping=False)
            x = ((S.cpu().numpy() + 1) // 2).astype(np.int64); f, phi = decode_eval(x, K, csr); pk.append(float(f)); ph.append(float(phi))
        sb.reset_env(); res["pkg-bSB"] = pk; res["pkg-phi"] = ph
        res["SA-native|64000"] = [float(sa_native(rng.integers(0, K, n), K, *csr, 2.0, 0.02, 64000, 300 + k)[1]) for k in range(4)]
        append_jsonl(out, dict(n=n, seed=seed, K=K, scale=float(W0[np.triu_indices(n, 1)].sum() / 3), res=res))
        print("Q", n, seed, {k: round(min(v), 3) for k, v in res.items()}, flush=True)

"""Ablation between the SB package's ballistic integrator and ours.
A numpy re-implementation of the package's step (simulated_bifurcation 2.0.0, symplectic_integrator.py and
simulated_bifurcation_optimizer.py), with one difference from ours switched at a time:
  init : oscillators uniform in [-1,1] (package) or [-0.1,0.1] (ours)
  order: momentum(pump) -> position -> momentum(coupling) -> walls (package) or momentum(pump+coupling) -> position -> walls (ours)
  aux  : field embedded as couplings to an extra oscillator (package) or applied as a constant field (ours)
  scale: package rule 0.5 sqrt(N)/||T||_F or ours 0.5 sqrt(N-1)/||couplings||_F
K = 3, 200 cells, seeds 0-4, eight agents, time step 0.5, 3,000 steps, linear pump. -> results/pkg_ablate.json"""
import json, math, itertools, warnings, numpy as np, torch, simulated_bifurcation as sb
warnings.filterwarnings("ignore")
from exp_common import *

K = 3


def run(T, steps, R, rng, init, order, aux, scale, ts=0.5):
    N = T.shape[0] - 1; C = T[:N, :N]; f = T[:N, N]
    c0 = 0.5 * math.sqrt(N) / math.sqrt((T ** 2).sum()) if scale == "pkg" else 0.5 * math.sqrt(N - 1) / math.sqrt((C ** 2).sum())
    w = 1.0 if init == "pkg" else 0.1; D = N + 1 if aux == "pkg" else N
    X = rng.uniform(-w, w, (D, R)); Y = rng.uniform(-w, w, (D, R))
    def force(X):            # package convention: minimise -0.5 x^T T x, so the force is +T x
        if aux == "pkg": return T @ X
        return C @ X + f[:, None]                     # auxiliary spin frozen at +1
    for t in range(steps):
        p = min(t / steps, 1.0)
        if order == "pkg":
            Y += ts * (p - 1.0) * X; X += ts * Y; Y += ts * c0 * force(X)
        else:
            Y += ts * ((p - 1.0) * X + c0 * force(X)); X += ts * Y
        m = np.abs(X) > 1; Y[m] = 0.0; np.clip(X, -1, 1, out=X)
    S = np.where(X >= 0, 1.0, -1.0)
    if aux == "pkg": S = S[-1][None, :] * S[:-1]
    return ((S.T + 1) // 2).astype(np.int64)


res = {}
for s in range(5):
    W, csr, a, J, lam, jj = setup(200, s, K)
    Jm = (J * 0.25).toarray(); h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
    T = sb.build_model(torch.tensor(Jm / 2, dtype=torch.float64), torch.tensor(h, dtype=torch.float64)).to_ising(domain="spin") \
          .as_simulated_bifurcation_tensor().numpy()
    for combo in itertools.product(("pkg", "ours"), repeat=4):
        X = run(T, 3000, 8, np.random.default_rng(s), *combo)
        r = [decode_eval(x, K, csr) for x in X]
        res.setdefault("|".join(combo), []).append((float(min(q for q, _ in r)), float(np.mean([p for _, p in r]))))
    print(s, flush=True)
summ = {k: dict(best_mean=float(np.mean([v[0] for v in L])), phi=float(np.mean([v[1] for v in L]))) for k, L in res.items()}
for k, v in sorted(summ.items(), key=lambda kv: kv[1]["phi"]): print("init|order|aux|scale =", k, {q: round(x, 3) for q, x in v.items()})
json.dump(dict(runs=res, summary=summ), open("results/pkg_ablate.json", "w"), indent=1)

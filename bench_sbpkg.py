"""
bench_sbpkg.py -- Experiment C: independent SB implementation (PyPI 'simulated-bifurcation' 2.0.0, PyTorch, CPU, 1 thread)
on the same Ising instances (rho = 0.25). Settings chosen on held-out seeds 100-101 (sbpkg_pilot.py), then frozen:
  ballistic        : time_step 0.1, pressure_slope 0.003, 3333 steps
  discrete, heated : time_step 0.25, pressure_slope 0.001, 4000 steps
128 agents, early stopping off, best agent by Ising energy, decoded with the same repair. 3 runs per instance.
"""
import sys, time, warnings, numpy as np, torch, simulated_bifurcation as sb
warnings.filterwarnings("ignore")
from exp_common import *
torch.set_num_threads(1)
CONF = {"pkg-bSB": dict(mode="ballistic", heated=False, ts=0.1, ps=0.003, steps=3333),
        "pkg-dSB-heated": dict(mode="discrete", heated=True, ts=0.25, ps=0.001, steps=4000)}

if __name__ == "__main__":
    K = int(sys.argv[1]); sizes = [int(x) for x in sys.argv[2].split(",")]; outf = sys.argv[3]
    done = done_keys(outf)
    for n in sizes:
        for seed in range(5):
            if (n, seed, K) in done: continue
            W, csr, a, J, lam, jj = setup(n, seed, K)
            Jm = (J * 0.25).toarray(); h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
            Q = torch.tensor(Jm / 2, dtype=torch.float32); l = torch.tensor(h, dtype=torch.float32)
            res = {}
            for name, c in CONF.items():
                sb.set_env(time_step=c["ts"], pressure_slope=c["ps"]); recs = []
                for k in range(3):
                    torch.manual_seed(1000 * seed + k); t0 = time.perf_counter()
                    S, E = sb.minimize(Q, l, domain="spin", agents=128, max_steps=c["steps"], mode=c["mode"], heated=c["heated"],
                                       verbose=False, best_only=True, early_stopping=False)
                    t = time.perf_counter() - t0; x = ((S.cpu().numpy() + 1) // 2).astype(np.int64)
                    f, phi = decode_eval(x, K, csr); recs.append(dict(f=float(f), phi=float(phi), t=t))
                res[name] = recs
            sb.reset_env()
            append_jsonl(outf, dict(n=n, seed=seed, K=K, runs=res))
            print(K, n, seed, {k: round(min(r["f"] for r in v), 2) for k, v in res.items()}, flush=True)

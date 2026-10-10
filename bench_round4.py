"""Round-4 checks on one machine.
A) Protocol P1 cell at 100 cells (K=3, 6 s per solver and instance) re-run with common seeds and recorded restart counts:
   native SA, and SA-QUBO / bSB at every penalty scale of Table XI.  -> results/p1_100.jsonl
B) Where the open-source SB package breaks: package (held-out settings and our settings) for n in {50, 100, 200} and
   agents in {8, 32, 128}, K=3, seeds 0-2: raw feasibility and Ising energy relative to our bSB.  -> results/sbpkg_bisect.jsonl"""
import sys, time, json, warnings, numpy as np
warnings.filterwarnings("ignore")
from exp_common import *
from bench_ran_qubo import sa_qubo


def part_a():
    out = "results/p1_100.jsonl"; done = done_keys(out, ("n", "seed", "K")); K = 3; n = 100; B = 6.0
    for seed in range(5):
        if (n, seed, K) in done: continue
        W = make_instance(n, seed, K); csr = csr_from_W(W); N = n * K; res = {}
        def loop(fn):
            rng = np.random.default_rng(7000 + seed); t0 = time.perf_counter(); best = 1e300; k = 0; fe = []
            while time.perf_counter() - t0 < B:
                e, f = fn(rng, k); best = min(best, e); fe.append(f); k += 1
            return dict(E=float(best), restarts=k, feas=float(np.mean(fe)))
        res["SA-native"] = loop(lambda rng, k: (sa_native(rng.integers(0, K, n), K, *csr, 2.0, 0.02, CFG["native_sweeps"], 70000 + k)[1], 1.0))
        for ls in (1.0, 0.5, 0.25, 0.12, 0.06):
            a, J, lam = build_qubo(W, K, ls); jp, jx, jd = J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64)
            res[f"SA-QUBO|{ls}"] = loop(lambda rng, k: decode_eval(sa_qubo(rng.integers(0, 2, N), a, jp, jx, jd, 0.5 * lam, 0.15, CFG["sweeps"], 71000 + k)[0], K, csr))
            def sbf(rng, k):
                Xb = sb_batch(a, J, CFG["sb_steps"], CFG["sb_R"], rng, False); es = [qubo_energy(z, a, jp, jx, jd) for z in Xb]
                return decode_eval(Xb[int(np.argmin(es))], K, csr)
            res[f"bSB|{ls}"] = loop(sbf)
        append_jsonl(out, dict(n=n, seed=seed, K=K, budget=B, res=res)); print("A", seed, {k: (round(v["E"], 2), v["restarts"]) for k, v in res.items()}, flush=True)


def part_b():
    import torch, simulated_bifurcation as sbp
    from bench_core import ising_of, sb_ising, ising_energy
    torch.set_num_threads(1); out = "results/sbpkg_bisect.jsonl"; K = 3
    for n in (50, 100, 200):
        for seed in range(3):
            W, csr, a, J, lam, jj = setup(n, seed, K); Jm, h = ising_of(a, J); rng = np.random.default_rng(9900 + seed)
            S = sb_ising(Jm, h, 3000, 8, rng, 0.5); E = ising_energy(S, Jm, h); e_ours = float(E.min())
            Q = torch.tensor(Jm.toarray() / 2, dtype=torch.float32); l = torch.tensor(h, dtype=torch.float32)
            for setting, (ts, ps, steps) in {"held-out": (0.1, 0.003, 3333), "ours": (0.5, 1 / 1500, 3000)}.items():
                for ag in (8, 32, 128):
                    sbp.set_env(time_step=ts, pressure_slope=ps); torch.manual_seed(100 * seed + ag)
                    s_, e_ = sbp.minimize(Q, l, domain="spin", agents=ag, max_steps=steps, mode="ballistic", verbose=False, best_only=True, early_stopping=False)
                    s_ = s_.numpy().astype(int); f, phi = decode_eval(((s_ + 1) // 2).astype(np.int64), K, csr)
                    rec = dict(n=n, seed=seed, K=K, setting=setting, agents=ag, phi=float(phi), E_pkg=float(ising_energy(s_[None, :], Jm, h)[0]),
                               E_ours=e_ours, E_allminus=float(ising_energy(-np.ones((1, len(h)), int), Jm, h)[0]), mean_spin=float(s_.mean()))
                    sbp.reset_env(); append_jsonl(out, rec)
            print("B", n, seed, flush=True)


if __name__ == "__main__":
    warm(3); part_a(); part_b(); print("DONE", flush=True)

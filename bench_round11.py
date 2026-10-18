"""(a) bSB at 256 times its baseline steps (768,000 steps), four runs per instance on the five P4 instances
of Table 6, to test whether the plateau continues; (b) protocol P1 for native SA (frozen schedule, restarts of 1,000 sweeps for 6 s)
repeated with ten seeds on the five 100-cell instances, to measure the seed-to-seed variation behind the 2-point reading threshold.
Usage: python3 bench_round11.py BSB256 K n a-b out.jsonl | SAREP n a-b reps out.jsonl"""
import sys, time, numpy as np
from exp_common import *

if __name__ == "__main__":
    mode = sys.argv[1]; warm(3); warm(4)
    if mode == "BSB256":
        K, n = int(sys.argv[2]), int(sys.argv[3]); a_, b_ = map(int, sys.argv[4].split("-")); outf = sys.argv[5]
        done = done_keys(outf)
        for seed in range(a_, b_):
            if (n, seed, K) in done: continue
            W, csr, a, J, lam, jj = setup(n, seed, K); rng = np.random.default_rng(110000 + 100 * K + seed); runs = []
            for k in range(4):
                t0 = time.perf_counter()
                Xb = sb_batch(a, J, 256 * CFG["sb_steps"], CFG["sb_R"], rng, False)
                es = [qubo_energy(x, a, *jj) for x in Xb]; x = Xb[int(np.argmin(es))]
                v, phi = decode_eval(x, K, csr)[:2]
                runs.append(dict(f=float(v), phi=float(phi), t=time.perf_counter() - t0))
            append_jsonl(outf, dict(n=n, seed=seed, K=K, runs=runs))
            print(K, n, seed, [round(r["f"], 2) for r in runs], [round(r["t"], 1) for r in runs], flush=True)
    elif mode == "SAREP":
        n = int(sys.argv[2]); a_, b_ = map(int, sys.argv[3].split("-")); reps = int(sys.argv[4]); outf = sys.argv[5]; K = 3
        done = done_keys(outf, ("n", "seed", "rep"))
        for rep in range(reps):
            for seed in range(a_, b_):
                if (n, seed, rep) in done: continue
                W, csr, a, J, lam, jj = setup(n, seed, K); rng = np.random.default_rng(111000 + 1000 * rep + seed)
                best = np.inf; cnt = 0; t0 = time.perf_counter()
                while time.perf_counter() - t0 < 6.0:
                    v = sa_native(rng.integers(0, K, n), K, *csr, 2.0, 0.02, CFG["native_sweeps"], int(rng.integers(1, 2**31 - 1)))[1]
                    best = min(best, float(v)); cnt += 1
                append_jsonl(outf, dict(n=n, seed=seed, rep=rep, K=K, f=best, restarts=cnt))
            print("rep", rep, "done", flush=True)

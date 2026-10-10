"""Native SA with its starting temperature tuned on held-out seeds 100-101 (bench_core.py part A; grid as for the polish):
T0 = 0.6 at K = 3 and 0.15 at K = 4 (lowest mean gap over 1k-64k sweeps), T1 = 0.02. Same effort plan, seeds and machine
as bench_hybrid.py, so it can be added to the matched-wall-clock analysis as 'SA-native (tuned)'."""
import sys, time, numpy as np
from exp_common import *
SA_PLAN = [(1, 40), (4, 20), (16, 10), (64, 5), (256, 3)]; T0 = {3: 0.6, 4: 0.15}
if __name__ == "__main__":
    K = int(sys.argv[1]); n = int(sys.argv[2]); outf = sys.argv[3]; a_, b_ = map(int, sys.argv[4].split("-"))
    warm(K); done = done_keys(outf)
    for seed in range(a_, b_):
        if (n, seed, K) in done: continue
        W, csr, a, J, lam, jj = setup(n, seed, K); rng = np.random.default_rng(43000 + seed); out = {}
        for f, nr in SA_PLAN:
            recs = []
            for k in range(nr):
                t0 = time.perf_counter(); v = sa_native(rng.integers(0, K, n), K, *csr, T0[K], 0.02, f * CFG["native_sweeps"], 910000 + 1000 * f + k)[1]
                recs.append(dict(f=float(v), t=time.perf_counter() - t0))
            out[str(f)] = recs
        append_jsonl(outf, dict(n=n, seed=seed, K=K, T0=T0[K], best=min(r["f"] for v in out.values() for r in v), runs=out))
        print(K, n, seed, flush=True)

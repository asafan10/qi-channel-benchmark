"""Reuse pattern + native-SA polish (adds to the reference): for every instance used anywhere in the paper, start native SA
from the best reuse pattern (reuse_baseline.py) and anneal at low temperature: T0 in {0.3, 0.15} -> 0.02, 64k sweeps, 2 runs each.
-> results/reuse_polish.jsonl (best cost per run)."""
import json, math, itertools, numpy as np
from exp_common import *


def best_reuse(n, K, csr):
    cols = int(math.ceil(math.sqrt(n * 1.15))); idx = np.arange(n); c = idx % cols; r = idx // cols; q = c - (r // 2); pats = []
    if K == 3:
        for a, b in itertools.product(range(3), range(3)):
            if (a, b) != (0, 0): pats.append((a * q + b * r) % 3)
    else:
        for a, b, a2, b2 in itertools.product(range(2), repeat=4):
            p = 2 * ((a * q + b * r) % 2) + ((a2 * q + b2 * r) % 2)
            if len(set(p.tolist())) == 4: pats.append(p)
    return min((p.astype(np.int64) for p in pats), key=lambda p: interference(p, *csr))


if __name__ == "__main__":
    out = "results/reuse_polish.jsonl"; done = done_keys(out)
    keys = [(K, n, s) for K, ns in ((3, (20, 50, 100, 200, 500)), (4, (50, 100, 200))) for n in ns for s in range(5)]
    keys += [(K, n, s) for K, n in ((3, 200), (3, 500), (4, 200)) for s in range(5, 30)]
    keys += [(K, 200, s) for K in (3, 4) for s in (100, 101, 102)]
    warm(3); warm(4)
    for K, n, s in keys:
        if (n, s, K) in done: continue
        W, csr, *_ = setup(n, s, K); p = best_reuse(n, K, csr); f0 = interference(p, *csr); runs = []
        for T0 in (0.3, 0.15):
            for k in range(2):
                ch, f = sa_native(p, K, *csr, T0, 0.02, 64000, 123000 + 10 * k + int(T0 * 100)); runs.append(dict(T0=T0, f=float(f)))
        append_jsonl(out, dict(K=K, n=n, seed=s, reuse=float(f0), best=float(min([f0] + [r["f"] for r in runs])), runs=runs))
        print(K, n, s, round(f0, 3), round(min(r["f"] for r in runs), 3), flush=True)

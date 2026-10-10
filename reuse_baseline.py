"""Assign channels by a regular reuse pattern on the hexagonal lattice indices,
ignoring the jitter. Cell i sits at column c = i mod cols and row r = i div cols (odd rows shifted right by half a cell);
axial coordinates q = c - floor(r/2), r. K = 3: the reuse-3 colourings (a q + b r) mod 3; K = 4: the reuse-4 pattern
2 (q mod 2) + (r mod 2) and its variants. The best pattern per instance is reported. -> results/reuse.json"""
import json, math, itertools, numpy as np
from exp_common import setup, interference
from refs import build_ref
ref = build_ref(); out = {}
for K, ns in ((3, (20, 50, 100, 200, 500)), (4, (50, 100, 200))):
    for n in ns:
        cols = int(math.ceil(math.sqrt(n * 1.15)))
        for s in range(5):
            W, csr, *_ = setup(n, s, K); best = None
            idx = np.arange(n); c = idx % cols; r = idx // cols; q = c - (r // 2)
            pats = []
            if K == 3:
                for a, b in itertools.product(range(3), range(3)):
                    if (a, b) != (0, 0): pats.append((a * q + b * r) % 3)
            else:
                for a, b, a2, b2 in itertools.product(range(2), repeat=4):
                    p = 2 * ((a * q + b * r) % 2) + ((a2 * q + b2 * r) % 2)
                    if len(set(p.tolist())) == 4: pats.append(p)
            for p in pats:
                f = interference(p.astype(np.int64), *csr)
                if best is None or f < best: best = f
            out[f"{K}_{n}_{s}"] = (best / ref[(K, n, s)] - 1) * 100
        print(K, n, round(float(np.mean([out[f"{K}_{n}_{s}"] for s in range(5)])), 1), flush=True)
json.dump(out, open("results/reuse.json", "w"), indent=1)

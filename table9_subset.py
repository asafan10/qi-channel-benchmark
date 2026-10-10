"""Table 9 restricted to seeds 0-4 (the instances of Table 7), and pool sizes of the expected best-of-N estimator in P6.
-> results/table9_subset.json"""
import json, io, contextlib, numpy as np
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    exec(open("analyse_hybrid.py").read())
out = {"seeds04": {}, "pools": {}}
for (K, n, s, T), (seeds, per) in PER.items():
    if s not in ("SA-native", "bSB"): continue
    idx = [i for i, x in enumerate(seeds) if x < 5]
    best = min((np.median(v[idx]) for e, v in per.items() if not np.isnan(v[idx]).any()), default=None)
    out["seeds04"][f"{K}|{n}|{s}|{T}"] = None if best is None else float(best)
for K, n in groups:
    sel = [r for r in rows if r["K"] == K and r["n"] == n]
    for s in ("SA-native", "SA-tuned", "bSB", "bSB+polish"):
        if s == "SA-tuned" and not all(r.get("sat") for r in sel): continue
        for e, (f, t) in samples(sel[0], s).items():
            pool = min(len(samples(r, s)[e][0]) for r in sel); tm = np.mean([samples(r, s)[e][1] for r in sel])
            out["pools"][f"{K}|{n}|{s}|{e}"] = dict(pool=int(pool), N_at_100s=float(100 // tm))
json.dump(out, open("results/table9_subset.json", "w"), indent=1)
for k, v in out["seeds04"].items():
    if k.startswith("3|200|SA-native") or k.startswith("3|200|bSB"): print(k, round(v, 2) if v is not None else v)
for k, v in out["pools"].items(): print(k, v)

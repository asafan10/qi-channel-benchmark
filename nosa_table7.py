"""Table 7 scored against the reference built without native SA's runs (refs.build_ref(no_native_sa=True)),
same pools, estimator and eligibility rule as analyse_fresh67.py -> results/nosa_table7.json"""
import json, glob, numpy as np
from scipy.special import gammaln
from refs import build_ref

RN = build_ref(no_native_sa=True); RF = build_ref(); BUD = [1, 3, 10, 30, 100]


def exp_min_exact(f, N):
    f = np.sort(np.asarray(f)); m = len(f)
    if N < 1 or 2 * N > m: return np.nan
    k = np.arange(1, m + 1); ok = (m - k) >= (N - 1); lw = np.full(m, -np.inf)
    lw[ok] = (gammaln(m - k[ok] + 1) - gammaln(N) - gammaln(m - k[ok] - N + 2)) - (gammaln(m + 1) - gammaln(N + 1) - gammaln(m - N + 1))
    return float(np.dot(np.exp(lw), f))


rows = {}
for f in glob.glob("results/fresh67_K*.jsonl"):
    for l in open(f): r = json.loads(l); rows[(r["K"], r["n"], r["seed"])] = r
for l in open("results/fresh67_topup.jsonl"):
    t = json.loads(l)
    for e, v in t["runs"].items(): rows[(t["K"], t["n"], t["seed"])]["runs"][t["kind"]][e] += v
out = {}
for K, n in ((3, 200), (3, 500), (4, 200)):
    seeds = sorted(s for (k, nn, s) in rows if (k, nn) == (K, n)); out[f"{K}_{n}"] = {}
    out[f"{K}_{n}"]["ref_shift"] = float(np.median([(RN[(K, n, s)] / RF[(K, n, s)] - 1) * 100 for s in seeds]))
    for sv in ("SA-tuned", "SA-frozen", "bSB"):
        res = {}
        for T in BUD:
            best = None
            for e, _ in rows[(K, n, seeds[0])]["runs"][sv].items():
                v = []
                for s in seeds:
                    p = rows[(K, n, s)]["runs"][sv][e]; N = int(T // np.mean([x["t"] for x in p]))
                    v.append((exp_min_exact([x["f"] for x in p], N) / RN[(K, n, s)] - 1) * 100)
                if np.isnan(v).any(): continue
                if best is None or np.median(v) < best: best = float(np.median(v))
            if best is not None: res[str(T)] = best
        out[f"{K}_{n}"][sv] = res
json.dump(out, open("results/nosa_table7.json", "w"), indent=1)
print(json.dumps(out, indent=1))

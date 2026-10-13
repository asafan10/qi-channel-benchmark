"""Summary of bench_round10.py (Numba constraint-preserving bSB): agreement with the NumPy version, time
per step against bSB, planted colourings, and matched wall-clock gaps with the estimator and pools of Table 7 -> results/round10_summary.json."""
import json, glob, numpy as np
from refs import build_ref
from scipy.special import gammaln


def exp_min_exact(f, N):        # as in analyse_fresh67.py (Table 7)
    f = np.sort(np.asarray(f)); m = len(f)
    if N < 1 or 2 * N > m: return np.nan
    k = np.arange(1, m + 1); ok = (m - k) >= (N - 1); lw = np.full(m, -np.inf)
    lw[ok] = (gammaln(m - k[ok] + 1) - gammaln(N) - gammaln(m - k[ok] - N + 2)) - (gammaln(m + 1) - gammaln(N + 1) - gammaln(m - N + 1))
    return float(np.dot(np.exp(lw), f))


REF = build_ref(); BUD = [1, 3, 10, 30, 100]; S = {}; med = lambda v: float(np.median(v))
NB = json.load(open("results/round10_NB.json"))["rows"]
S["nb"] = {f"{K}_{st}": dict(gap=med([med([(f / r["ref"] - 1) * 100 for f in r["runs"][st]["f"]]) for r in NB if r["K"] == K]),
                             t=float(np.mean([r["runs"][st]["t"] for r in NB if r["K"] == K]))) for K in (3, 4) for st in ("3000", "12000")}
PL = json.load(open("results/round10_PL.json"))["rows"]
S["planted"] = {}
for n in (50, 100, 200):
    rr = [r for r in PL if r["n"] == n]
    for st in rr[0]["res"]:
        S["planted"][f"{n}_{st}"] = dict(hits=int(sum(min(r["res"][st]) == 0.0 for r in rr)), best_pct=med([min(r["res"][st]) / r["scale"] * 100 for r in rr]),
                                         median_pct=med([med(r["res"][st]) / r["scale"] * 100 for r in rr]))
cp = {}
for f in glob.glob("results/round10_CPPN_K*.jsonl"):
    for l in open(f): r = json.loads(l); cp[(r["K"], r["seed"])] = r
fr = {}
for f in glob.glob("results/fresh67_K*.jsonl"):
    for l in open(f): r = json.loads(l); fr[(r["K"], r["n"], r["seed"])] = r
for l in open("results/fresh67_topup.jsonl"):
    t = json.loads(l)
    for e, v in t["runs"].items(): fr[(t["K"], t["n"], t["seed"])]["runs"][t["kind"]][e] += v
S["cpp"] = {}
for K in (3, 4):
    seeds = sorted(s for (k, s) in cp if k == K); ref = {s: REF[(K, 200, s)] for s in seeds}
    pools = {"CP-bSB": {s: {e: [(x["f"], x["t"]) for x in v] for e, v in cp[(K, s)]["runs"].items()} for s in seeds}}
    for sv in ("SA-tuned", "bSB"):
        pools[sv] = {s: {e: [(x["f"], x["t"]) for x in v] for e, v in fr[(K, 200, s)]["runs"][sv].items()} for s in seeds}
    single = {e: dict(gap=med([med([(f / ref[s] - 1) * 100 for f, _ in pools["CP-bSB"][s][e]]) for s in seeds]),
                      t=float(np.mean([np.mean([t for _, t in pools["CP-bSB"][s][e]]) for s in seeds]))) for e in pools["CP-bSB"][seeds[0]]}
    bsb_t = float(np.mean([np.mean([t for _, t in pools["bSB"][s]["1"]]) for s in seeds]))
    per = {}
    for sv, P_ in pools.items():
        per[sv] = {}
        for T in BUD:
            best = None
            for e in P_[seeds[0]]:
                v = np.array([(exp_min_exact([f for f, _ in P_[s][e]], int(T // np.mean([t for _, t in P_[s][e]]))) / ref[s] - 1) * 100 for s in seeds])
                if np.isnan(v).any(): continue
                if best is None or np.median(v) < best[1]: best = (e, float(np.median(v)), v)
            if best: per[sv][T] = dict(effort=best[0], gap=best[1], v=best[2].tolist())
    common = [T for T in BUD if T in per["SA-tuned"] and T in per["CP-bSB"]]
    S["cpp"][str(K)] = dict(n=len(seeds), single=single, slowdown=single["3000"]["t"] / bsb_t,
                            matched={sv: {str(T): dict(effort=d["effort"], gap=d["gap"]) for T, d in per[sv].items()} for sv in per},
                            wins_SA={str(T): int(sum(a < b for a, b in zip(per["SA-tuned"][T]["v"], per["CP-bSB"][T]["v"]))) for T in common},
                            vs_bSB={str(T): int(sum(a < b for a, b in zip(per["CP-bSB"][T]["v"], per["bSB"][T]["v"]))) for T in common if T in per["bSB"]})
json.dump(S, open("results/round10_summary.json", "w"), indent=1)
print(json.dumps(S, indent=1)[:5000])

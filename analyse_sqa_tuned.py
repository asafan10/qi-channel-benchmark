"""Tuned SQA at matched wall-clock time (bench_sqa_tuned.py) against tuned native SA and bSB on the pools of Table 7, with the estimator and
eligibility rule of Table 7. -> results/sqa_tuned_summary.json (table rows)"""
import json, glob, numpy as np
from scipy.special import gammaln
from refs import build_ref

REF = build_ref(); BUD = [1, 3, 10, 30, 100]; med = lambda v: float(np.median(v))


def exp_min_exact(f, N):
    f = np.sort(np.asarray(f)); m = len(f)
    if N < 1 or 2 * N > m: return np.nan
    k = np.arange(1, m + 1); ok = (m - k) >= (N - 1); lw = np.full(m, -np.inf)
    lw[ok] = (gammaln(m - k[ok] + 1) - gammaln(N) - gammaln(m - k[ok] - N + 2)) - (gammaln(m + 1) - gammaln(N + 1) - gammaln(m - N + 1))
    return float(np.dot(np.exp(lw), f))


fr = {}
for f in glob.glob("results/fresh67_K*.jsonl"):
    for l in open(f): r = json.loads(l); fr[(r["K"], r["n"], r["seed"])] = r
for l in open("results/fresh67_topup.jsonl"):
    t = json.loads(l)
    for e, v in t["runs"].items(): fr[(t["K"], t["n"], t["seed"])]["runs"][t["kind"]][e] += v
TU = json.load(open("results/sqa_tune.json")); S = dict(tune={K: dict(best=TU[K]["best"], best_gap=TU[K]["best_gap"], default_gap=TU[K]["default_gap"]) for K in TU})
rows = []; f1 = lambda x: "–" if x is None else f"{x:.1f}".replace("-", "−")
for K in (3, 4):
    sq = {json.loads(l)["seed"]: json.loads(l) for l in open(f"results/sqa_pool_K{K}.jsonl")}; seeds = sorted(sq); ref = {s: REF[(K, 200, s)] for s in seeds}
    pools = {"SQA": {s: {e: [(x["f"], x["t"]) for x in v] for e, v in sq[s]["runs"].items()} for s in seeds}}
    for sv in ("SA-tuned", "bSB"):
        pools[sv] = {s: {e: [(x["f"], x["t"]) for x in v] for e, v in fr[(K, 200, s)]["runs"][sv].items()} for s in seeds}
    per = {}
    for sv, P_ in pools.items():
        per[sv] = {}
        for T in BUD:
            best = None
            for e in P_[seeds[0]]:
                v = np.array([(exp_min_exact([f for f, _ in P_[s][e]], int(T // np.mean([t for _, t in P_[s][e]]))) / ref[s] - 1) * 100 for s in seeds])
                if np.isnan(v).any(): continue
                if best is None or np.median(v) < best[1]: best = (e, float(np.median(v)), v)
            if best: per[sv][T] = best
    single = {e: med([med([(f / ref[s] - 1) * 100 for f, _ in pools["SQA"][s][e]]) for s in seeds]) for e in pools["SQA"][seeds[0]]}
    phi = min(x["phi"] for s in seeds for v in sq[s]["runs"].values() for x in v)
    wins = {T: int(sum(a < b for a, b in zip(per["SA-tuned"][T][2], per["SQA"][T][2]))) for T in BUD if T in per["SA-tuned"] and T in per["SQA"]}
    S[str(K)] = dict(n=len(seeds), single=single, phi_min=phi, wins_SA=wins, matched={sv: {str(T): d[1] for T, d in per[sv].items()} for sv in per})
    for T in BUD:
        rows.append([f"{K}", f"{T} s"] + [f1(per[sv][T][1]) if T in per[sv] else "–" for sv in ("SQA", "SA-tuned", "bSB")])
S["header"] = ["$K$", "Budget", "SQA (tuned)", "SA-native (tuned)", "bSB"]; S["rows"] = rows; S["widths"] = [600, 900, 1500, 1500, 1500]
c3, c4 = S["3"], S["4"]; t3, t4 = S["tune"]["3"], S["tune"]["4"]
m = lambda c, sv, T: c["matched"][sv].get(str(T))
json.dump(S, open("results/sqa_tuned_summary.json", "w"), indent=1, ensure_ascii=False)
[print(r) for r in rows]

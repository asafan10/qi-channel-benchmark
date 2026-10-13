"""Summary of bench_round9.py -> results/round9_summary.json. Matched wall-clock comparisons of the constraint-preserving bSB use the same estimator and eligibility rule as
Table 7 (exact expected best of N without replacement, pool >= 2N runs on every instance) and the fresh pools of Table 7 on the same machine."""
import json, glob, os, numpy as np
from scipy.special import gammaln
from refs import build_ref

REF = build_ref(); BUD = [1, 3, 10, 30, 100]; rng = np.random.default_rng(69); S = {}
med = lambda v: float(np.median(v))


def exp_min_exact(f, N):        # as in analyse_fresh67.py
    f = np.sort(np.asarray(f)); m = len(f)
    if N < 1 or 2 * N > m: return np.nan
    k = np.arange(1, m + 1); ok = (m - k) >= (N - 1); lw = np.full(m, -np.inf)
    lw[ok] = (gammaln(m - k[ok] + 1) - gammaln(N) - gammaln(m - k[ok] - N + 2)) - (gammaln(m + 1) - gammaln(N + 1) - gammaln(m - N + 1))
    return float(np.dot(np.exp(lw), f))


# ---- REP: decoding choices
R = json.load(open("results/round9_REP.json"))["rows"]; S["rep"] = {}
for K in (3, 4):
    for st in ("3000", "12000", "48000"):
        rr = [r for r in R if r["K"] == K]
        S["rep"][f"{K}_{st}"] = {k: med([med([(b[k] / r["ref"] - 1) * 100 for b in r["runs"][st]]) for r in rr]) for k in ("default", "all8", "visited_E", "visited_best")}
        S["rep"][f"{K}_{st}"]["same_state"] = float(np.mean([b["same_state"] for r in rr for b in r["runs"][st]]))
R4 = json.load(open("results/round9_REP4.json"))["rows"]; S["rep4"] = {}
for lab in ("row", "perm"):
    rr = [r for r in R4 if r["labels"] == lab]
    S["rep4"][lab] = dict(default=med([(x["default"] / r["ref"] - 1) * 100 for r in rr for x in r["runs"]]),
                          visited=med([(x["vis_f"] / r["ref"] - 1) * 100 for r in rr for x in r["runs"]]),
                          visited_range=[min((x["vis_f"] / r["ref"] - 1) * 100 for r in rr for x in r["runs"]), max((x["vis_f"] / r["ref"] - 1) * 100 for r in rr for x in r["runs"])],
                          vis_phi_max=max(x["vis_phi"] for r in rr for x in r["runs"]), vis_phi_min=min(x["vis_phi"] for r in rr for x in r["runs"]),
                          vis_step=[min(x["vis_step"] for r in rr for x in r["runs"]), max(x["vis_step"] for r in rr for x in r["runs"])])
# ---- OH2: one-hot bSB at K = 2 against the Max-Cut encoding of round 6
O = json.load(open("results/round9_OH2.json"))["rows"]; K2 = json.load(open("results/round6_summary.json"))["K2"]["0.15"]
S["oh2"] = {st: dict(gap=med([med([(x["f"] / r["best_round6"] - 1) * 100 for x in r["runs"][st]["rows"]]) for r in O]),
                     phi_min=min(x["phi"] for r in O for x in r["runs"][st]["rows"]), t=float(np.mean([r["runs"][st]["t"] for r in O]))) for st in ("3000", "12000", "48000")}
S["oh2_maxcut_bsb"] = {f: K2["bsb"][f]["median_gap"] for f in ("1", "4", "16")}; S["oh2_sa"] = {f: K2["sa"][f]["median_gap"] for f in ("4", "16", "64", "256")}
# ---- POLT
P = json.load(open("results/round9_POLT.json"))["rows"]
S["polt"] = {k: dict(t=med([p[k]["t"] for p in P]), acc=med([p[k]["acc"] for p in P])) for k in ("frozen_from_random", "polish_from_bsb")}
# ---- CPT / CPTB / CPP: constraint-preserving bSB
G = json.load(open("results/round9_CPT.json")); GB = json.load(open("results/round9_CPTB.json"))
S["cpt"] = {K: dict(best=min((v for v in G.values() if v["K"] == K), key=lambda v: v["mean_gap"]), default_bsb=GB[f"K{K}"]["mean_gap"],
                    grid=[min(v["mean_gap"] for v in G.values() if v["K"] == K), max(v["mean_gap"] for v in G.values() if v["K"] == K)]) for K in (3, 4)}
cp = {}
for f in glob.glob("results/round9_CPP_K*.jsonl"):
    for l in open(f): r = json.loads(l); cp[(r["K"], r["seed"])] = r
fr = {}
for f in glob.glob("results/fresh67_K*.jsonl"):
    for l in open(f): r = json.loads(l); fr[(r["K"], r["n"], r["seed"])] = r
if os.path.exists("results/fresh67_topup.jsonl"):
    for l in open("results/fresh67_topup.jsonl"):
        t = json.loads(l)
        for e, v in t["runs"].items(): fr[(t["K"], t["n"], t["seed"])]["runs"][t["kind"]][e] += v
S["cpp"] = {}
for K in (3, 4):
    seeds = sorted(s for (k, s) in cp if k == K)
    if not seeds: continue
    ref = {s: REF[(K, 200, s)] for s in seeds}
    pools = {"CP-bSB": {s: {e: [(x["f"], x["t"]) for x in v] for e, v in cp[(K, s)]["runs"].items()} for s in seeds}}
    for sv in ("SA-tuned", "bSB"):
        pools[sv] = {s: {e: [(x["f"], x["t"]) for x in v] for e, v in fr[(K, 200, s)]["runs"][sv].items()} for s in seeds}
    single = {e: dict(gap=med([med([(f / ref[s] - 1) * 100 for f, _ in pools["CP-bSB"][s][e]]) for s in seeds]),
                      t=float(np.mean([np.mean([t for _, t in pools["CP-bSB"][s][e]]) for s in seeds]))) for e in pools["CP-bSB"][seeds[0]]}
    per = {}
    for sv, P_ in pools.items():
        per[sv] = {}
        for T in BUD:
            best = None
            for e in P_[seeds[0]]:
                v = []
                for s in seeds:
                    p = P_[s][e]; tt = np.mean([t for _, t in p]); N = int(T // tt)
                    v.append((exp_min_exact([f for f, _ in p], N) / ref[s] - 1) * 100)
                v = np.array(v)
                if np.isnan(v).any(): continue
                if best is None or np.median(v) < best[1]: best = (e, float(np.median(v)), v)
            if best: per[sv][T] = dict(effort=best[0], gap=best[1], v=best[2].tolist())
    wins = {str(T): int(sum(a < b for a, b in zip(per["SA-tuned"][T]["v"], per["CP-bSB"][T]["v"]))) for T in BUD if T in per["SA-tuned"] and T in per["CP-bSB"]}
    S["cpp"][str(K)] = dict(n=len(seeds), single=single, matched={sv: {str(T): {k: d[k] for k in ("effort", "gap")} for T, d in per[sv].items()} for sv in per},
                            wins_SA=wins)
json.dump(S, open("results/round9_summary.json", "w"), indent=1)
print(json.dumps({k: S[k] for k in S if k not in ("rep",)}, indent=1)[:6000])

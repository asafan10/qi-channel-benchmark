"""Summary of bench_round7.py (C0, POL, K4T, HEAT, HEATP) -> results/round7_summary.json.
The heated-bSB pools (HEATP) are compared at matched wall-clock time with the pools of Table 7 (bench_fresh67.py, same machine, same
instances: K = 3, 200 cells, seeds 0-4) using the same exact expected best of N and the same eligibility rule (pool >= 2N runs)."""
import json, glob, numpy as np
from refs import build_ref
from scipy.special import gammaln
rng = np.random.default_rng(67)


def exp_min_exact(f, N):        # as in analyse_fresh67.py
    f = np.sort(np.asarray(f)); m = len(f)
    if N < 1 or 2 * N > m: return np.nan
    k = np.arange(1, m + 1); ok = (m - k) >= (N - 1); lw = np.full(m, -np.inf)
    lw[ok] = (gammaln(m - k[ok] + 1) - gammaln(N) - gammaln(m - k[ok] - N + 2)) - (gammaln(m + 1) - gammaln(N + 1) - gammaln(m - N + 1))
    return float(np.dot(np.exp(lw), f))


def boot_med(v, B=4000):
    v = np.asarray(v); b = [np.median(rng.choice(v, len(v), replace=True)) for _ in range(B)]
    return float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))

REF = build_ref(); BUD = [1, 3, 10, 30, 100]; out = {}
H = json.load(open("results/round7_HEAT.json"))
out["heat"] = dict(chosen=H["heat"], tune=H["tune"],
                   single={k: float(np.median([(np.median(r["runs"][k]) / r["ref"] - 1) * 100 for r in H["rows"]])) for k in H["rows"][0]["runs"]})
K4 = json.load(open("results/round7_K4T_eval.json"))
out["k4t"] = dict(best=K4["best"], cfg=K4["cfg"], phi_checked=K4["phi_checked"], default_grid=K4["default_grid"], grid_best=K4["grid_best"],
                  single={k: float(np.median([(np.median(r["runs"][k]) / r["ref"] - 1) * 100 for r in K4["rows"]])) for k in K4["rows"][0]["runs"]})
rows = []
for f in sorted(glob.glob("results/round7_HEATP_*.json")): rows += json.load(open(f))["rows"]
fr = {}
for f in glob.glob("results/fresh67_K3_200.jsonl"):
    for l in open(f): r = json.loads(l); fr[r["seed"]] = r
seeds = sorted(r["seed"] for r in rows); HP = {r["seed"]: r for r in rows}


def pools(sv, s):
    if sv == "bSB-heated": return {e: [(x["f"], x["t"]) for x in v] for e, v in HP[s]["runs"].items()}
    return {e: [(x["f"], x["t"]) for x in v] for e, v in fr[s]["runs"][sv].items()}


per = {}
for sv in ("SA-tuned", "bSB", "bSB-heated"):
    per[sv] = {}
    for T in BUD:
        best = None
        for e in pools(sv, seeds[0]):
            v = []
            for s in seeds:
                p = pools(sv, s)[e]; tt = np.mean([t for _, t in p]); N = int(T // tt)
                v.append((exp_min_exact([f for f, _ in p], N) / REF[(3, 200, s)] - 1) * 100)
            v = np.array(v)
            if np.isnan(v).any(): continue
            if best is None or np.median(v) < best[1]: best = (e, float(np.median(v)), v)
        if best: lo, hi = boot_med(best[2]); per[sv][T] = dict(effort=best[0], gap=best[1], lo=lo, hi=hi, v=best[2].tolist())
out["matched"] = {sv: {str(T): {k: v for k, v in d.items() if k != "v"} for T, d in per[sv].items()} for sv in per}
out["wins_SA_vs_heated"] = {str(T): int(sum(a < b for a, b in zip(per["SA-tuned"][T]["v"], per["bSB-heated"][T]["v"])))
                            for T in BUD if T in per["SA-tuned"] and T in per["bSB-heated"]}
out["heated_single"] = {e: dict(gap=float(np.median([np.median([(x["f"] / REF[(3, 200, s)] - 1) * 100 for x in HP[s]["runs"][e]]) for s in seeds])),
                             t=float(np.mean([np.mean([x["t"] for x in HP[s]["runs"][e]]) for s in seeds])),
                             phi=float(np.mean([np.mean([x["phi"] for x in HP[s]["runs"][e]]) for s in seeds]))) for e in HP[seeds[0]]["runs"]}
out["n_seeds"] = len(seeds)
json.dump(out, open("results/round7_summary.json", "w"), indent=1)
print(json.dumps({k: out[k] for k in ("matched", "wins_SA_vs_heated", "heated_single", "k4t")}, indent=1)[:3000])

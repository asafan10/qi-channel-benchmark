"""Analysis of bench_sigma.py (protocol P6 on the sigma = 0.3 layout, 200 cells, K = 3), with the conventions of analyse_hybrid.py:
expected best-of-N at matched wall-clock budgets (best effort per budget, median over instances, bootstrap CI), paired Wilcoxon
tests with Holm correction and median paired differences with bootstrap CIs, single-run plateau of bSB, gaps to the SDP bound.
Reference per instance: best plan of any run on it (including reuse and reuse + SA). Output results/sigma03.json."""
import json, numpy as np
from scipy.stats import wilcoxon

MARGIN = 0.5   # equivalence margin in gap points (two one-sided Wilcoxon tests)
BUDGETS = [1, 3, 10, 30, 100]; rng = np.random.default_rng(3)
R = [json.loads(l) for l in open("results/sigma03.jsonl")]; R.sort(key=lambda r: r["seed"])


def samples(r, which):
    d = {}
    if which in ("SA-native", "SA-tuned"):
        src = r["hybrid"]["SA"] if which == "SA-native" else r["tuned"]
        for f, recs in src.items(): d[f"x{f}"] = (np.array([x["f"] for x in recs]), np.mean([x["t"] for x in recs]))
        return d
    for f, recs in r["hybrid"]["bSB"].items():
        if which == "bSB": d[f"x{f}"] = (np.array([x["f_sel"] for x in recs]), np.mean([x["t"] for x in recs]))
        else:
            for s in recs[0]["pol"]: d[f"x{f}+{int(s)//1000}k"] = (np.array([x["pol"][s]["f"] for x in recs]), np.mean([x["t"] + x["pol"][s]["t"] for x in recs]))
    if which == "bSB+polish":
        for s in r["ext"][0]["pol"]: d[f"x1+{int(s)//1000}k (ext)"] = (np.array([x["pol"][s]["f"] for x in r["ext"]]), np.mean([x["t"] + x["pol"][s]["t"] for x in r["ext"]]))
    return d


def allcosts(r):
    out = [r["reuse"], r["reuse_sa"]]
    for w in ("SA-native", "SA-tuned", "bSB", "bSB+polish"):
        for f, _ in samples(r, w).values(): out += list(f)
    return min(out)


import os
ref = {r["seed"]: allcosts(r) for r in R}; lb = {r["seed"]: r["lb"] for r in R}
SDPREF = os.environ.get("SDPREF") == "1"     # robustness: reference also includes the SDP-rounding plans (sdp_round.py)
if SDPREF:
    for l in open("results/sdp_round_sigma0.3.jsonl"):
        d_ = json.loads(l); ref[d_["seed"]] = min(ref[d_["seed"]], d_["f_best"])
OUTF = "results/sigma03_sdpref_v1.json" if SDPREF else "results/sigma03_v1.json"   # first-round analysis; the manuscript uses analyse_sigma_ext.py


def exp_min(f, N, B=300):
    N = min(N, 1000); return float(np.mean([rng.choice(f, size=N, replace=True).min() for _ in range(B)]))


SOL = ["SA-native", "SA-tuned", "bSB", "bSB+polish"]; res = {}; INST = {}; PER = {}
for s in SOL:
    for T in BUDGETS:
        per = {}
        for r in R:
            for e, (f, t) in samples(r, s).items():
                N = int(T // t); per.setdefault(e, []).append((exp_min(f, N) / ref[r["seed"]] - 1) * 100 if N >= 1 else np.nan)
        PER[(s, T)] = ({e: np.array(v) for e, v in per.items()})
        best = None
        for e, v in per.items():
            v = np.array(v)
            if np.isnan(v).any(): continue
            if best is None or np.median(v) < best[1]: best = (e, float(np.median(v)), v)
        if best is None: continue
        bs = [np.median(rng.choice(best[2], size=len(best[2]))) for _ in range(2000)]
        res[f"{s}|{T}"] = [best[1], *np.percentile(bs, [2.5, 97.5]).tolist(), best[0]]; INST[(s, T)] = best[2]


def tests(pairs):
    out = []
    for a_, b_, T in pairs:
        if (a_, T) not in INST or (b_, T) not in INST: continue
        dl = INST[(a_, T)] - INST[(b_, T)]; r2 = np.random.default_rng(T)
        bs = [np.median(r2.choice(dl, size=len(dl))) for _ in range(2000)]
        p_t = max(wilcoxon(dl + MARGIN, alternative="greater").pvalue, wilcoxon(dl - MARGIN, alternative="less").pvalue)
        out.append(dict(a=a_, b=b_, T=T, p=float(wilcoxon(dl).pvalue) if np.any(dl != 0) else 1.0, med_diff=float(np.median(dl)), p_tost=float(p_t),
                        ci90=np.percentile(bs, [5, 95]).tolist(),
                        ci=np.percentile(bs, [2.5, 97.5]).tolist(), wins_a=int((dl < 0).sum()), N=len(dl)))
    order = sorted(range(len(out)), key=lambda i: out[i]["p"]); mx = 0.0
    for k, i in enumerate(order): mx = max(mx, min(1.0, (len(out) - k) * out[i]["p"])); out[i]["p_holm"] = mx
    order = sorted(range(len(out)), key=lambda i: out[i]["p_tost"]); mx = 0.0
    for k, i in enumerate(order): mx = max(mx, min(1.0, (len(out) - k) * out[i]["p_tost"])); out[i]["p_tost_holm"] = mx
    return out


t_pol = tests([(a, "bSB+polish", T) for T in BUDGETS for a in ("SA-native", "SA-tuned")])
t_sb = tests([("SA-native", "bSB", T) for T in BUDGETS])
single = {f: float(np.median([np.median((np.array([x["f_sel"] for x in r["hybrid"]["bSB"][f]]) / ref[r["seed"]] - 1) * 100) for r in R]))
          for f in R[0]["hybrid"]["bSB"]}
sa_single = {f: float(np.median([np.median((np.array([x["f"] for x in r["hybrid"]["SA"][f]]) / ref[r["seed"]] - 1) * 100) for r in R]))
             for f in R[0]["hybrid"]["SA"]}
out = dict(n_inst=len(R), res=res, tests_vs_polish=t_pol, tests_sa_vs_bsb=t_sb, bsb_single=single, sa_single=sa_single,
           reuse_gap=float(np.mean([(r["reuse"] / ref[r["seed"]] - 1) * 100 for r in R])),
           cert_gap=float(np.mean([(ref[s] / lb[s] - 1) * 100 for s in ref])), cert_gap_max=float(max((ref[s] / lb[s] - 1) * 100 for s in ref)),
           sdp_gap_100={s: float(np.median([(ref[r["seed"]] * (1 + g / 100)) / lb[r["seed"]] * 100 - 100 for r, g in zip(R, INST[(s, 100)])])) for s in SOL if (s, 100) in INST})
seeds = np.array([r["seed"] for r in R])
def crossfit(key):
    per = {e: v for e, v in PER[key].items() if not np.isnan(v).any()}; A = seeds < 15; o = np.full(len(seeds), np.nan)
    for tr, te in ((A, ~A), (~A, A)):
        e_ = min(per, key=lambda e: np.median(per[e][tr])); o[te] = per[e_][te]
    return o
var = {"crossfit_gap": [], "crossfit_raw": [], "crossfit_sdp": []}
for a_, b_ in (("SA-native", "bSB+polish"), ("SA-tuned", "bSB+polish"), ("SA-native", "bSB")):
    for T in BUDGETS:
        ga, gb = crossfit((a_, T)), crossfit((b_, T))
        ra = np.array([ref[x] * (1 + g / 100) for x, g in zip(seeds, ga)]); rb = np.array([ref[x] * (1 + g / 100) for x, g in zip(seeds, gb)])
        for nm, dl in (("crossfit_gap", ga - gb), ("crossfit_raw", ra - rb), ("crossfit_sdp", (ra - rb) / np.array([lb[x] for x in seeds]) * 100)):
            var[nm].append(dict(a=a_, b=b_, T=T, p=float(wilcoxon(dl).pvalue), med=float(np.median(dl)),
                                p_tost=float(max(wilcoxon(dl + MARGIN, alternative="greater").pvalue, wilcoxon(dl - MARGIN, alternative="less").pvalue)) if nm != "crossfit_raw" else None))
for nm, tt in var.items():
    for fam in (lambda t: t["b"] == "bSB+polish", lambda t: t["b"] == "bSB"):
        sub = [t for t in tt if fam(t)]; o = sorted(range(len(sub)), key=lambda i: sub[i]["p"]); mx = 0.0
        for k, i in enumerate(o): mx = max(mx, min(1.0, (len(sub) - k) * sub[i]["p"])); sub[i]["p_holm"] = mx
        if nm != "crossfit_raw":
            o = sorted(range(len(sub)), key=lambda i: sub[i]["p_tost"]); mx = 0.0
            for k, i in enumerate(o): mx = max(mx, min(1.0, (len(sub) - k) * sub[i]["p_tost"])); sub[i]["p_tost_holm"] = mx
out["variants"] = var
for nm, tt in var.items():
    for t in tt: print(nm, t["a"], t["b"], t["T"], "med %.3f holm %.3g" % (t["med"], t["p_holm"]), "" if t.get("p_tost_holm") is None else "tost %.3g" % t["p_tost_holm"])
json.dump(out, open(OUTF, "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ("tests_vs_polish", "tests_sa_vs_bsb")}, indent=0))
for t in t_pol + t_sb: print(t["a"], "vs", t["b"], t["T"], "med %.2f [%.2f, %.2f]" % (t["med_diff"], *t["ci"]), "90%% [%.2f, %.2f]" % tuple(t["ci90"]), "holm %.3g" % t["p_holm"], "tost %.3g holm %.3g" % (t["p_tost"], t["p_tost_holm"]), "wins %d/%d" % (t["wins_a"], t["N"]))

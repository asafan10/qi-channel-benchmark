"""Re-analysis of protocol P6 on the irregular layout with the fresh pools of bench_sigma_ext.py, following PREREG_sigma03_ext.md.
Exact expected best-of-N without replacement, efforts eligible only where m >= 2N on every instance; paired Wilcoxon and TOST
(margin 0.5 points) with Holm over the five budgets. Also: equal-time comparison of semidefinite rounding with tuned SA (exploratory).
-> results/sigma03_ext.json"""
import json, sys
import numpy as np
from math import comb
from scipy.stats import wilcoxon

MARGIN = 0.5; BUDGETS = [1, 3, 10, 30, 100]; BOOT = 2000
OLD = {r["seed"]: r for r in map(json.loads, open("results/sigma03.jsonl"))}
EXT = {r["seed"]: r for r in map(json.loads, open("results/sigma03_ext.jsonl"))}
import os as _os
if _os.path.exists("results/sigma03_ext2.jsonl") and _os.environ.get("NO_EXT2") != "1":   # amendment 3: extended pools
    for r2 in map(json.loads, open("results/sigma03_ext2.jsonl")):
        R0 = EXT[r2["seed"]]["runs"]
        for f, recs in r2["runs"]["bSB"].items(): R0["bSB"][f] = R0["bSB"][f] + recs
        for key in ("SA-tuned", "SA-frozen"):
            for f, recs in r2["runs"][key].items(): R0[key][f] = R0[key][f] + recs
SDP = {r["seed"]: r for r in map(json.loads, open("results/sdp_round_sigma0.3.jsonl"))}
import os
SEEDS = sorted(s for s in EXT if s in OLD)
if os.environ.get("EXCL_EARLY") == "1":                 # post hoc sensitivity: drop the four instances run before the pause, whose timings were perturbed
    SEEDS = [s for s in SEEDS if s not in (0, 1, 15, 16)]


def old_costs(r):
    v = [r["reuse"], r["reuse_sa"]]
    for recs in r["hybrid"]["SA"].values(): v += [x["f"] for x in recs]
    for recs in r["tuned"].values(): v += [x["f"] for x in recs]
    for recs in r["hybrid"]["bSB"].values():
        for x in recs: v += [x["f_sel"], x.get("f_desc", np.inf)] + [p["f"] for p in x["pol"].values()]
    for x in r["ext"]: v += [p["f"] for p in x["pol"].values()]
    return v


def samples(r, which):
    R = r["runs"]; d = {}
    if which in ("SA-tuned", "SA-frozen"):
        for f, recs in R[which].items(): d[f"x{f}"] = [(x["f"], x["t"]) for x in recs]
    elif which == "bSB":
        for f, recs in R["bSB"].items(): d[f"x{f}"] = [(x["f_sel"], x["t"]) for x in recs]
    else:
        for s in R["bSB"]["1"][0]["pol"]:
            d[f"x1+{int(s) // 1000}k"] = [(x["pol"][s]["f"], x["t"] + x["pol"][s]["t"]) for x in R["bSB"]["1"]]
    return d


def new_costs(r):
    v = []
    for w_ in ("SA-tuned", "SA-frozen", "bSB", "bSB+polish"):
        for runs in samples(r, w_).values(): v += [f for f, _ in runs]
    return v


REF = {s: min(old_costs(OLD[s]) + new_costs(EXT[s])) for s in SEEDS}
REF_SDP = {s: min(REF[s], SDP[s]["f_best"]) if s in SDP else REF[s] for s in SEEDS}
LB = {s: OLD[s]["lb"] for s in SEEDS}


def exp_min_exact(f, N):
    f = np.sort(np.asarray(f)); m = len(f)
    if N < 1 or N > m: return np.nan
    w = np.array([comb(m - k, N - 1) for k in range(1, m + 1)], dtype=float) / comb(m, N)
    return float(np.dot(w, f))


SOL = ["SA-tuned", "SA-frozen", "bSB", "bSB+polish"]
def table(ref):
    res, per = {}, {}
    for s_ in SOL:
        effs = samples(EXT[SEEDS[0]], s_).keys()
        for T in BUDGETS:
            cand = {}
            for e in effs:
                vals = []; ok = True
                for sd in SEEDS:
                    runs = samples(EXT[sd], s_)[e]; t = np.mean([x[1] for x in runs]); N = int(T // t); m = len(runs)
                    if N < 1 or m < 2 * N: ok = False; break
                    vals.append((exp_min_exact([x[0] for x in runs], N) / ref[sd] - 1) * 100)
                if ok: cand[e] = np.array(vals)
            if not cand: res[f"{s_}|{T}"] = None; continue
            e_ = min(cand, key=lambda e: np.median(cand[e])); v = cand[e_]
            rng = np.random.default_rng(T); bs = [np.median(rng.choice(v, len(v))) for _ in range(BOOT)]
            res[f"{s_}|{T}"] = [float(np.median(v)), *np.percentile(bs, [2.5, 97.5]).tolist(), e_]; per[(s_, T)] = v
    return res, per


def holm(ps):
    o = sorted(range(len(ps)), key=lambda i: ps[i]); out = [0.0] * len(ps); mx = 0.0
    for k, i in enumerate(o): mx = max(mx, min(1.0, (len(ps) - k) * ps[i])); out[i] = mx
    return out


def tests(per, a_, b_, ref):
    rows = []
    for T in BUDGETS:
        if (a_, T) not in per or (b_, T) not in per: continue
        dl = per[(a_, T)] - per[(b_, T)]
        rng = np.random.default_rng(100 + T); bs = [np.median(rng.choice(dl, len(dl))) for _ in range(BOOT)]
        p = float(wilcoxon(dl).pvalue) if np.any(dl != 0) else 1.0
        pt = float(max(wilcoxon(dl + MARGIN, alternative="greater").pvalue, wilcoxon(dl - MARGIN, alternative="less").pvalue))
        rows.append(dict(a=a_, b=b_, T=T, med_diff=float(np.median(dl)), ci95=np.percentile(bs, [2.5, 97.5]).tolist(),
                         ci90=np.percentile(bs, [5, 95]).tolist(), wins_a=int((dl < 0).sum()), N=len(dl), p=p, p_tost=pt,
                         ratio_med=float(np.median(per[(a_, T)]) / np.median(per[(b_, T)])) if np.median(per[(b_, T)]) > 0 else None))
    for k, h in zip(rows, holm([r["p"] for r in rows])): k["p_holm"] = h
    for k, h in zip(rows, holm([r["p_tost"] for r in rows])): k["p_tost_holm"] = h
    return rows


OUT = {"seeds": SEEDS, "pools": {}, "eligibility": {}, "timing_by_process": {}}
for s_ in SOL:
    for e in samples(EXT[SEEDS[0]], s_):
        for T in BUDGETS:
            cnt = 0; Ns = []
            for sd in SEEDS:
                runs = samples(EXT[sd], s_)[e]; N = int(T // np.mean([x[1] for x in runs])); Ns.append(N)
                cnt += int(N >= 1 and len(runs) >= 2 * N)
            OUT["eligibility"][f"{s_}|{e}|{T}"] = dict(eligible_instances=cnt, N_min=int(min(Ns)), N_max=int(max(Ns)),
                                                        m=min(len(samples(EXT[sd], s_)[e]) for sd in SEEDS))
        ta = [np.mean([x[1] for x in samples(EXT[sd], s_)[e]]) for sd in SEEDS if sd < 15]; tb = [np.mean([x[1] for x in samples(EXT[sd], s_)[e]]) for sd in SEEDS if sd >= 15]
        OUT["timing_by_process"][f"{s_}|{e}"] = dict(seeds_0_14=float(np.mean(ta)) if ta else None, seeds_15_29=float(np.mean(tb)) if tb else None)
for s_ in SOL:
    for e, runs in samples(EXT[SEEDS[0]], s_).items():
        OUT["pools"][f"{s_}|{e}"] = dict(m=min(len(samples(EXT[sd], s_)[e]) for sd in SEEDS),
                                         t=float(np.mean([np.mean([x[1] for x in samples(EXT[sd], s_)[e]]) for sd in SEEDS])))
for name, ref in (("primary", REF), ("sdp_ref", REF_SDP)):
    res, per = table(ref)
    T_pol = tests(per, "SA-tuned", "bSB+polish", ref); T_fro = tests(per, "SA-frozen", "bSB+polish", ref); T_bsb = tests(per, "SA-tuned", "bSB", ref)
    eq = [t["T"] for t in T_pol if t["p_tost_holm"] < 0.05]; tested = [t["T"] for t in T_pol]; upto = 0
    for T in BUDGETS:
        if T in tested and T in eq: upto = T      # amendment 2: stop at the first untested or non-equivalent budget
        else: break
    OUT[name] = dict(res=res, tuned_vs_polish=T_pol, frozen_vs_polish=T_fro, tuned_vs_bsb=T_bsb, equivalent_up_to=upto,
                     cert_gap=float(np.mean([(ref[s] / LB[s] - 1) * 100 for s in SEEDS])))
    if name == "primary": PER = per

# single-run plateau of bSB alone (fresh pools)
OUT["bsb_single"] = {e: float(np.median([np.median([(f / REF[sd] - 1) * 100 for f, _ in samples(EXT[sd], "bSB")[e]]) for sd in SEEDS]))
                     for e in samples(EXT[SEEDS[0]], "bSB")}
OUT["bsb_best_of_pool"] = {e: float(np.median([(min(f for f, _ in samples(EXT[sd], "bSB")[e]) / REF[sd] - 1) * 100 for sd in SEEDS]))
                           for e in samples(EXT[SEEDS[0]], "bSB")}

# semidefinite rounding: ablation, variability and an equal-time comparison with tuned SA (exploratory)
sd_ = [s for s in SEEDS if s in SDP]
g = lambda v, s: (v / REF[s] - 1) * 100
abl = dict(descent_best=[g(SDP[s]["f_desc"], s) for s in sd_], polish_best=[g(SDP[s]["f_pol"], s) for s in sd_],
           descent_median=[g(SDP[s]["f_desc_med"], s) for s in sd_], polish_median=[g(1 + SDP[s]["polish_median"] / 100, s) for s in sd_],
           t_desc=[SDP[s]["t_sdp"] + SDP[s]["t_round"] for s in sd_], t_total=[SDP[s]["t_sdp"] + SDP[s]["t_round"] + SDP[s]["t_polish"] for s in sd_],
           eps=sorted({SDP[s]["eps"] for s in sd_}))
OUT["rounding"] = {k: (v if k == "eps" else dict(median=float(np.median(v)), min=float(np.min(v)), max=float(np.max(v)))) for k, v in abl.items()}
eq_rows = []
for label, fkey, tkey in (("rounding + descent", "descent_best", "t_desc"), ("rounding + descent + polish", "polish_best", "t_total")):
    dl, ok = [], True
    for i, s in enumerate(sd_):
        T = abl[tkey][i]; best = None
        for e, runs in samples(EXT[s], "SA-tuned").items():
            t = np.mean([x[1] for x in runs]); N = int(T // t); m = len(runs)
            if N < 1 or m < 2 * N: continue
            v = g(exp_min_exact([x[0] for x in runs], N), s)
            best = v if best is None else min(best, v)          # most favourable effort for SA on each instance
        dl.append(abl[fkey][i] - best)
    dl = np.array(dl)
    eq_rows.append(dict(method=label, med_diff=float(np.median(dl)), wins_rounding=int((dl < 0).sum()), N=len(dl), p=float(wilcoxon(dl).pvalue)))
OUT["rounding_equal_time"] = eq_rows
json.dump(OUT, open("results/sigma03_ext_excl_early.json" if os.environ.get("EXCL_EARLY") == "1" else "results/sigma03_ext.json", "w"), indent=1)
for name in ("primary", "sdp_ref"):
    print("==", name, "equivalent up to", OUT[name]["equivalent_up_to"], "s; cert gap %.2f" % OUT[name]["cert_gap"])
    for k, v in OUT[name]["res"].items(): print("  ", k, None if v is None else "%.2f [%.2f, %.2f] %s" % tuple(v))
    for t in OUT[name]["tuned_vs_polish"] + OUT[name]["frozen_vs_polish"] + OUT[name]["tuned_vs_bsb"]:
        print("  ", t["a"], "vs", t["b"], t["T"], "med %.2f ci90 [%.2f, %.2f] wins %d/%d holm %.3g tostH %.3g ratio %s" % (
            t["med_diff"], *t["ci90"], t["wins_a"], t["N"], t["p_holm"], t["p_tost_holm"], t["ratio_med"]))
print("bSB single", OUT["bsb_single"], "best of pool", OUT["bsb_best_of_pool"])
print("rounding", json.dumps(OUT["rounding"])); print("equal time", OUT["rounding_equal_time"])


# ---------------------------------------------------------------- outputs in the schema of analyse_sigma.py, so the manuscript uses the fresh pools
NAMEMAP = {"SA-frozen": "SA-native", "SA-tuned": "SA-tuned", "bSB": "bSB", "bSB+polish": "bSB+polish"}
def eligible(s_, T):
    out = {}
    for e in samples(EXT[SEEDS[0]], s_):
        ok = True; Ns = []
        for sd in SEEDS:
            runs = samples(EXT[sd], s_)[e]; N = int(T // np.mean([x[1] for x in runs]))
            if N < 1 or len(runs) < 2 * N: ok = False; break
            Ns.append(N)
        if ok: out[e] = Ns
    return out


def gaps(s_, e, T, ref):
    return np.array([(exp_min_exact([x[0] for x in samples(EXT[sd], s_)[e]], int(T // np.mean([x[1] for x in samples(EXT[sd], s_)[e]]))) / ref[sd] - 1) * 100 for sd in SEEDS])


def crossfit(s_, T, ref):
    el = eligible(s_, T); seeds = np.array(SEEDS); A = seeds < 15; o = np.full(len(seeds), np.nan)
    per = {e: gaps(s_, e, T, ref) for e in el}
    for tr, te in ((A, ~A), (~A, A)):
        e_ = min(per, key=lambda e: np.median(per[e][tr])); o[te] = per[e_][te]
    return o


def schema(ref, res, per):
    o = {"n_inst": len(SEEDS), "res": {}, "estimator": "exact expected best-of-N without replacement, efforts with m >= 2N"}
    for k, v in res.items():
        s_, T = k.split("|"); o["res"][f"{NAMEMAP[s_]}|{T}"] = v
    def conv(rows):
        out = []
        for r in rows:
            r = dict(r); r["a"] = NAMEMAP[r["a"]]; r["b"] = NAMEMAP[r["b"]]; r["ci"] = r.pop("ci95"); out.append(r)
        return out
    o["tests_vs_polish"] = conv(tests(per, "SA-frozen", "bSB+polish", ref)) + conv(tests(per, "SA-tuned", "bSB+polish", ref))
    o["tests_sa_vs_bsb"] = conv(tests(per, "SA-frozen", "bSB", ref))
    o["bsb_single"] = {e[1:]: float(np.median([np.median([(f / ref[sd] - 1) * 100 for f, _ in samples(EXT[sd], "bSB")[e]]) for sd in SEEDS]))
                       for e in samples(EXT[SEEDS[0]], "bSB")}
    o["sa_single"] = {e[1:]: float(np.median([np.median([(f / ref[sd] - 1) * 100 for f, _ in samples(EXT[sd], "SA-frozen")[e]]) for sd in SEEDS]))
                      for e in samples(EXT[SEEDS[0]], "SA-frozen")}
    o["reuse_gap"] = float(np.mean([(OLD[s]["reuse"] / ref[s] - 1) * 100 for s in SEEDS]))
    o["cert_gap"] = float(np.mean([(ref[s] / LB[s] - 1) * 100 for s in SEEDS])); o["cert_gap_max"] = float(max((ref[s] / LB[s] - 1) * 100 for s in SEEDS))
    o["sdp_gap_100"] = {NAMEMAP[s_]: float(np.median([ref[sd] * (1 + g_ / 100) / LB[sd] * 100 - 100 for sd, g_ in zip(SEEDS, per[(s_, 100)])]))
                        for s_ in SOL if (s_, 100) in per}
    o["sdp_gap"] = {f"{NAMEMAP[s_]}|{T}": float(np.median([ref[sd] * (1 + g_ / 100) / LB[sd] * 100 - 100 for sd, g_ in zip(SEEDS, per[(s_, T)])]))
                    for s_ in SOL for T in BUDGETS if (s_, T) in per}
    var = {"crossfit_gap": [], "crossfit_raw": [], "crossfit_sdp": []}
    for a_, b_ in (("SA-frozen", "bSB+polish"), ("SA-tuned", "bSB+polish"), ("SA-frozen", "bSB")):
        for T in BUDGETS:
            if not eligible(a_, T) or not eligible(b_, T): continue
            ga, gb = crossfit(a_, T, ref), crossfit(b_, T, ref)
            ra = np.array([ref[x] * (1 + g_ / 100) for x, g_ in zip(SEEDS, ga)]); rb = np.array([ref[x] * (1 + g_ / 100) for x, g_ in zip(SEEDS, gb)])
            for nm, dl in (("crossfit_gap", ga - gb), ("crossfit_raw", ra - rb), ("crossfit_sdp", (ra - rb) / np.array([LB[x] for x in SEEDS]) * 100)):
                var[nm].append(dict(a=NAMEMAP[a_], b=NAMEMAP[b_], T=T, p=float(wilcoxon(dl).pvalue), med=float(np.median(dl)),
                                    p_tost=float(max(wilcoxon(dl + MARGIN, alternative="greater").pvalue, wilcoxon(dl - MARGIN, alternative="less").pvalue)) if nm != "crossfit_raw" else None))
    for nm, tt in var.items():
        for key in set((t["a"], t["b"]) for t in tt):
            sub = [t for t in tt if (t["a"], t["b"]) == key]
            for t, h in zip(sub, holm([t["p"] for t in sub])): t["p_holm"] = h
            if nm != "crossfit_raw":
                for t, h in zip(sub, holm([t["p_tost"] for t in sub])): t["p_tost_holm"] = h
    o["variants"] = var
    # descriptive: smallest margin (0.01 grid) at which equivalence holds at each tested budget, Holm over the tested budgets
    curve = {}
    dls = {T: per[("SA-tuned", T)] - per[("bSB+polish", T)] for T in BUDGETS if ("SA-tuned", T) in per and ("bSB+polish", T) in per}
    for k in range(1, 101):
        m_ = k / 100; Ts = sorted(dls); ps = [float(max(wilcoxon(dls[T] + m_, alternative="greater").pvalue, wilcoxon(dls[T] - m_, alternative="less").pvalue)) for T in Ts]
        for T, h in zip(Ts, holm(ps)):
            if h < 0.05 and str(T) not in curve: curve[str(T)] = m_
    o["margin_curve"] = curve
    v1 = json.load(open("results/sigma03_v1.json"))
    o["first_round"] = {str(t["T"]): t["med_diff"] for t in v1["tests_vs_polish"] if t["a"] == "SA-tuned"}
    tb = [v for v in OUT["timing_by_process"].values() if v["seeds_0_14"] and v["seeds_15_29"]]
    dr = max(abs(v["seeds_0_14"] / v["seeds_15_29"] - 1) for v in tb) * 100 if tb else float("nan")
    o["timing_drift_pct"] = dr
    return o


for name, ref, fn in ((() if os.environ.get("EXCL_EARLY") == "1" else (("primary", REF, "results/sigma03.json"), ("sdp_ref", REF_SDP, "results/sigma03_sdpref.json")))):
    res, per = table(ref)
    json.dump(schema(ref, res, per), open(fn, "w"), indent=1)

# semidefinite rounding on the irregular layout, against the fresh reference, for Section 5.10 and Table B9
SDJ = json.load(open("results/sdp_round.json")) if os.environ.get("EXCL_EARLY") != "1" else None
gs = np.array([g(SDP[s]["f_best"], s) for s in sd_])
if SDJ is not None: SDJ["sigma"] = dict(N=len(sd_), gap_ref_median=float(np.median(gs)), gap_ref_min=float(gs.min()), gap_ref_max=float(gs.max()),
                    n_improved=int((gs < -1e-7).sum()), gap_lb_mean=float(np.mean([(SDP[s]["f_best"] / LB[s] - 1) * 100 for s in sd_])),
                    ref_lb_mean=float(np.mean([(REF[s] / LB[s] - 1) * 100 for s in sd_])), t_median=OUT["rounding"]["t_total"]["median"],
                    t_max=OUT["rounding"]["t_total"]["max"], eps=OUT["rounding"]["eps"], ablation=OUT["rounding"], equal_time=OUT["rounding_equal_time"])
if SDJ is not None: json.dump(SDJ, open("results/sdp_round.json", "w"), indent=1)

"""Analysis of Experiment A (bench_hybrid.py): single-run quality and matched wall-clock (Eq. 25) for
SA-native, bSB, bSB+descent and bSB+SA-polish. Reference = min over every run of every experiment available for the instance
(this file, plus final_reference.json / earlier jsonl files if present)."""
import json, sys, glob, numpy as np
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BUDGETS = [1, 3, 10, 30, 100]; BOOT = 2000; rng = np.random.default_rng(3)
rows = []
for f in sys.argv[1:] or sorted(glob.glob("results/hybrid_K*.jsonl")):
    rows += [json.loads(l) for l in open(f)]
from refs import build_ref
ref = build_ref()
ext = {}
for f in sorted(glob.glob("results/hybext_K*.jsonl")):
    for l in open(f):
        e = json.loads(l); k = (e["K"], e["n"], e["seed"]); ext[k] = e["runs"]
for r in rows: r["ext"] = ext.get((r["K"], r["n"], r["seed"]))
sat = {}
for f in sorted(glob.glob("results/satuned_K*.jsonl")):
    for l in open(f):
        e = json.loads(l); sat[(e["K"], e["n"], e["seed"])] = e["runs"]
for r in rows: r["sat"] = sat.get((r["K"], r["n"], r["seed"]))
for f in glob.glob("results/ref_extra_*.json"):           # optional extra references {"K_n_seed": cost}
    for k, v in json.load(open(f)).items():
        K, n, s = map(int, k.split("_")); ref[(K, n, s)] = min(ref.get((K, n, s), np.inf), v)


def samples(r, which):
    """-> {effort label: (array of costs, mean time per run)}"""
    d = {}
    if which == "SA-tuned":
        for f, recs in (r.get("sat") or {}).items():
            d[f"x{f}"] = (np.array([x["f"] for x in recs]), np.mean([x["t"] for x in recs]))
    elif which == "SA-native":
        for f, recs in r["runs"]["SA"].items():
            d[f"x{f}"] = (np.array([x["f"] for x in recs]), np.mean([x["t"] for x in recs]))
    else:
        for f, recs in r["runs"]["bSB"].items():
            if which == "bSB":
                d[f"x{f}"] = (np.array([x["f_sel"] for x in recs]), np.mean([x["t"] for x in recs]))
            elif which == "bSB+descent":
                d[f"x{f}"] = (np.array([x["f_desc"] for x in recs]), np.mean([x["t"] + x["t_desc"] for x in recs]))
            elif which == "bSB+polish":
                for s in recs[0]["pol"]:
                    d[f"x{f}+{int(s)//1000}k"] = (np.array([x["pol"][s]["f"] for x in recs]), np.mean([x["t"] + x["pol"][s]["t"] for x in recs]))
        if which == "bSB+polish" and r.get("ext"):
            for s_ in r["ext"][0]["pol"]:
                d[f"x1+{int(s_)//1000}k (ext)"] = (np.array([x["pol"][s_]["f"] for x in r["ext"]]), np.mean([x["t"] + x["pol"][s_]["t"] for x in r["ext"]]))
    return d


def exp_min(f, N, B=300):
    N = min(N, 1000); return float(np.mean([rng.choice(f, size=N, replace=True).min() for _ in range(B)]))


SOLVERS = ["SA-native", "SA-tuned", "bSB", "bSB+descent", "bSB+polish"]
PER = {}; INST = {}   # (K, n, solver, T) -> (effort, seeds, per-instance values)
groups = sorted({(r["K"], r["n"]) for r in rows})
lines = ["# Experiment A: native-space refinement after bSB decoding\n",
         "Reference: minimum over every run of every experiment available (refs.py). Same machine for all timings.\n"]
res = {}
for K, n in groups:
    sel = sorted([r for r in rows if r["K"] == K and r["n"] == n], key=lambda r: r["seed"])
    lines += [f"\n## K={K}, n={n} ({len(sel)} instances)\n", "### Single-run median gap (%) and mean time per run (ms)\n",
              "| solver | effort | median gap % | time/run (ms) | runs/instance |", "|---|---|---|---|---|"]
    for s in SOLVERS:
        if s == "SA-tuned" and not all(r.get("sat") for r in sel): continue
        for e in samples(sel[0], s):
            g = np.median([np.median((samples(r, s)[e][0] / ref[(K, n, r["seed"])] - 1) * 100) for r in sel])
            t = np.mean([samples(r, s)[e][1] for r in sel]); nr = len(samples(sel[0], s)[e][0])
            lines.append(f"| {s} | {e} | {g:.1f} | {1000*t:.0f} | {nr} |")
    lines += ["\n### Matched wall-clock: expected best-of-N gap (%), median over instances [95% bootstrap CI], best effort per budget\n",
              "| T | " + " | ".join(SOLVERS) + " |", "|---|" + "---|" * len(SOLVERS)]
    for T in BUDGETS:
        row = []
        for s in SOLVERS:
            per = {}
            if s == "SA-tuned" and not all(r.get("sat") for r in sel):
                res[(K, n, s, T)] = None; row.append("-"); continue
            for r in sel:
                for e, (f, t) in samples(r, s).items():
                    N = int(T // t)
                    per.setdefault(e, []).append((exp_min(f, N) / ref[(K, n, r["seed"])] - 1) * 100 if N >= 1 else np.nan)
            PER[(K, n, s, T)] = ([r["seed"] for r in sel], {e: np.array(v) for e, v in per.items()})
            best = None
            for e, v in per.items():
                v = np.array(v)
                if np.isnan(v).any(): continue
                m = np.median(v)
                if best is None or m < best[1]: best = (e, m, v)
            if best is None:
                res[(K, n, s, T)] = None; row.append("-"); continue
            bs = [np.median(rng.choice(best[2], size=len(best[2]), replace=True)) for _ in range(BOOT)]
            lo, hi = np.percentile(bs, [2.5, 97.5]); res[(K, n, s, T)] = (best[1], lo, hi, best[0]); INST[(K, n, s, T)] = (best[0], [r["seed"] for r in sel], best[2])
            row.append(f"{best[1]:.1f} [{lo:.1f}, {hi:.1f}] ({best[0]})")
        lines.append(f"| {T} s | " + " | ".join(row) + " |")
txt = "\n".join(lines) + "\n"; Path("results/hybrid.md").write_text(txt); print(txt)

# figure: same visual conventions as Figures 1-2 of the paper; 600 dpi; full text width
from figstyle import FONT as _F0
plt.rcParams.update(_F0)
fig, axs = plt.subplots(1, len(groups), figsize=(7.16, 2.95), squeeze=False)
from figstyle import STY as _FS, XLAB, YLAB, FONT, title
sty = {s: _FS["SA-frozen" if s == "SA-native" else s] for s in ("SA-native", "SA-tuned", "bSB", "bSB+descent", "bSB+polish")}
for ax, (K, n) in zip(axs[0], groups):
    for s in SOLVERS:
        pts = [(T, res[(K, n, s, T)]) for T in BUDGETS if res.get((K, n, s, T))]
        if not pts: continue
        x = [p[0] for p in pts]; m = [p[1][0] for p in pts]; lo = [p[1][1] for p in pts]; hi = [p[1][2] for p in pts]
        c, mk, ls, lab = sty[s]; ax.plot(x, m, color=c, marker=mk, ls=ls, label=lab, ms=4, lw=1.2); ax.fill_between(x, lo, hi, color=c, alpha=.15)
    ax.axvline(1, color="k", lw=1.0, ls="--", alpha=.5)
    ax.axhline(5, color="grey", lw=0.9, ls=":"); ax.text(95, 5.25, "5% gap", fontsize=7, ha="right", color="grey")
    ax.set_xscale("log"); ax.set_xlabel(XLAB)
    if ax is axs[0][0]: ax.set_ylabel(YLAB)
    ax.set_ylim(bottom=-0.5); ax.set_title(title(K, n)); ax.grid(alpha=.3, which="both")
h, l = axs[0][0].get_legend_handles_labels(); fig.legend(h, l, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.85), w_pad=0.8); fig.savefig("figures/hybrid_matched_wallclock.pdf"); fig.savefig("figures/hybrid_matched_wallclock.png", dpi=600)
json.dump({f"{K}|{n}|{s}|{T}": v for (K, n, s, T), v in res.items()}, open("results/hybrid_res.json", "w"))


# ---------------------------------------------------------------- robustness: paired tests, reference without native SA, SDP gaps
from scipy.stats import wilcoxon
from refs import sdp_bounds
ref_nosa = build_ref(no_native_sa=True); LB = sdp_bounds()
tests = []
for K, n in groups:
    for T in BUDGETS:
        for a_ in ("SA-native", "SA-tuned"):
            A = INST.get((K, n, a_, T)); P = INST.get((K, n, "bSB+polish", T))
            if not A or not P: continue
            dlt = np.array(A[2]) - np.array(P[2])
            p = wilcoxon(dlt).pvalue if np.any(dlt != 0) else 1.0
            _r2 = np.random.default_rng(1000 * K + n + T); _bs = [np.median(_r2.choice(dlt, size=len(dlt), replace=True)) for _ in range(2000)]
            _A = np.array([A[2][i] for i in range(len(dlt))]); _P = np.array(P[2])
            _pt = max(wilcoxon(dlt + 0.5, alternative="greater").pvalue, wilcoxon(dlt - 0.5, alternative="less").pvalue)
            tests.append(dict(K=K, n=n, T=T, a=a_, p=float(p), med_diff=float(np.median(dlt)), ci=[float(x) for x in np.percentile(_bs, [2.5, 97.5])],
                              ci90=[float(x) for x in np.percentile(_bs, [5, 95])], p_tost=float(_pt),
                              mean_gap_a=float(np.mean(_A)), mean_gap_pol=float(np.mean(_P)), raw_a=float(np.mean([ref[(K, n, sd)] * (1 + g / 100) for sd, g in zip(A[1], A[2])])), raw_pol=float(np.mean([ref[(K, n, sd)] * (1 + g / 100) for sd, g in zip(P[1], P[2])])), wins_a=int((dlt < 0).sum()), N=len(dlt)))
order = sorted(range(len(tests)), key=lambda i: tests[i]["p"]); m = len(tests); run_max = 0.0
for rank, i in enumerate(order):                     # Holm step-down
    run_max = max(run_max, min(1.0, (m - rank) * tests[i]["p"])); tests[i]["p_holm"] = run_max
order = sorted(range(len(tests)), key=lambda i: tests[i]["p_tost"]); run_max = 0.0
for rank, i in enumerate(order): run_max = max(run_max, min(1.0, (m - rank) * tests[i]["p_tost"])); tests[i]["p_tost_holm"] = run_max


def cost_at(K, n, s, T, r):
    eff = INST[(K, n, s, T)][0]; f, t = samples(r, s)[eff]; N = int(T // t)
    return exp_min(f, N) if N >= 1 else np.nan


nosa = {}; sdpg = {}
for K, n in groups:
    sel = sorted([r for r in rows if r["K"] == K and r["n"] == n], key=lambda r: r["seed"])
    for T in (10, 100):
        for s in ("SA-native", "SA-tuned", "bSB+polish"):
            if (K, n, s, T) not in INST: continue
            c = [cost_at(K, n, s, T, r) for r in sel]
            nosa[f"{K}|{n}|{s}|{T}"] = float(np.median([(ci / ref_nosa[(K, n, r["seed"])] - 1) * 100 for ci, r in zip(c, sel)]))
            lbs = [(ci, LB[(K, n, r["seed"])]) for ci, r in zip(c, sel) if (K, n, r["seed"]) in LB]
            if len(lbs) == 5: sdpg[f"{K}|{n}|{s}|{T}"] = float(np.median([(ci / l - 1) * 100 for ci, l in lbs]))
json.dump(dict(tests=tests, nosa=nosa, sdp_gap=sdpg, n_inst={f"{K}|{n}": len([r for r in rows if r["K"] == K and r["n"] == n]) for K, n in groups}),
          open("results/robust.json", "w"), indent=1)
print("robust: %d tests; nosa %s; sdp %s" % (len(tests), nosa, sdpg))


# ---------------------------------------------------------------- robustness of the paired tests: cross-fitted effort choice,
# raw-cost differences, and differences normalised by the semidefinite bound (200 cells, where all 30 bounds exist)
def crossfit(key):
    seeds, per = PER[key]; seeds = np.array(seeds); A = seeds < 15; out = np.full(len(seeds), np.nan)
    for train, test in ((A, ~A), (~A, A)):
        ok = {e: v for e, v in per.items() if not np.isnan(v).any()}
        e_star = min(ok, key=lambda e: np.median(ok[e][train])); out[test] = ok[e_star][test]
    return seeds, out
LB2 = dict(LB)
if Path("results/sdp_rig_p6.jsonl").exists():
    for l_ in open("results/sdp_rig_p6.jsonl"):
        q_ = json.loads(l_); LB2[(q_["K"], q_["n"], q_["seed"])] = q_["lb_f_rig"]
variants = {"crossfit_gap": [], "crossfit_raw": [], "crossfit_sdp": []}
for K, n in groups:
    for T in BUDGETS:
        for a_ in ("SA-native", "SA-tuned"):
            if (K, n, a_, T) not in PER or (K, n, "bSB+polish", T) not in PER: continue
            sd, ga = crossfit((K, n, a_, T)); _, gp = crossfit((K, n, "bSB+polish", T))
            raw_a = np.array([ref[(K, n, x)] * (1 + g / 100) for x, g in zip(sd, ga)]); raw_p = np.array([ref[(K, n, x)] * (1 + g / 100) for x, g in zip(sd, gp)])
            for name, dl in (("crossfit_gap", ga - gp), ("crossfit_raw", raw_a - raw_p)):
                variants[name].append(dict(K=K, n=n, T=T, a=a_, p=float(wilcoxon(dl).pvalue) if np.any(dl != 0) else 1.0, med=float(np.median(dl)), N=len(dl)))
            if all((K, n, x) in LB2 for x in sd):
                dl = np.array([(ra - rp) / LB2[(K, n, x)] * 100 for ra, rp, x in zip(raw_a, raw_p, sd)])
                variants["crossfit_sdp"].append(dict(K=K, n=n, T=T, a=a_, p=float(wilcoxon(dl).pvalue) if np.any(dl != 0) else 1.0, med=float(np.median(dl)), N=len(dl)))
for name, tt in variants.items():
    order = sorted(range(len(tt)), key=lambda i: tt[i]["p"]); mx = 0.0
    for rank, i in enumerate(order): mx = max(mx, min(1.0, (len(tt) - rank) * tt[i]["p"])); tt[i]["p_holm"] = mx
R_ = json.load(open("results/robust.json")); R_["variants"] = variants; json.dump(R_, open("results/robust.json", "w"), indent=1)
base = {(t["K"], t["n"], t["T"], t["a"]): t for t in R_["tests"]}
def verdict(t, sign): return "SA" if t["p_holm"] < 0.05 and sign < 0 else ("polish" if t["p_holm"] < 0.05 else "ns")
for name, tt in variants.items():
    changed = [(t["K"], t["n"], t["T"], t["a"], verdict(base[(t["K"], t["n"], t["T"], t["a"])], base[(t["K"], t["n"], t["T"], t["a"])]["med_diff"]), verdict(t, t["med"])) for t in tt
               if verdict(t, t["med"]) != verdict(base[(t["K"], t["n"], t["T"], t["a"])], base[(t["K"], t["n"], t["T"], t["a"])]["med_diff"])]
    print(name, len(tt), "tests; verdicts changed:", changed)

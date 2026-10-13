"""Tables 6 and 7 and Figures 1-2 from the fresh one-machine pools of bench_fresh67.py.
Table 6: median single-run gap (%) to the solver-only best-known cost and mean time per run, by effort.
Table 7: exact expected best of N = floor(T / t_e) runs without replacement, using only efforts whose pool holds at least 2N runs on
every instance; the effort with the lowest median over instances is reported with a 95% bootstrap interval over instances.
Paired Wilcoxon tests (tuned SA against untuned bSB, and tuned against frozen SA) at each budget, Holm within each setting.
-> results/fresh67_summary.json, figures/fig_K3.png, figures/fig_K4.png (and .pdf)"""
import json, glob, numpy as np
from scipy.special import gammaln
from scipy.stats import wilcoxon
from refs import build_ref
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

REF = build_ref(); BUD = [1, 3, 10, 30, 100]; rng = np.random.default_rng(67); BOOT = 4000
SOLVERS = ["SA-tuned", "SA-frozen", "bSB", "bSB-tuned"]
rows = {}
for f in glob.glob("results/fresh67_K*.jsonl"):
    for l in open(f):
        r = json.loads(l); rows[(r["K"], r["n"], r["seed"])] = r
import os
if os.path.exists("results/fresh67_topup.jsonl"):      # extra K = 4 SA runs (topup_fresh67.py)
    for l in open("results/fresh67_topup.jsonl"):
        t = json.loads(l)
        for e, v in t["runs"].items(): rows[(t["K"], t["n"], t["seed"])]["runs"][t["kind"]][e] += v
SETTINGS = [(3, 200), (3, 500), (4, 200)]


def exp_min_exact(f, N):
    f = np.sort(np.asarray(f)); m = len(f)
    if N < 1 or 2 * N > m: return np.nan
    k = np.arange(1, m + 1); ok = (m - k) >= (N - 1)
    lw = np.full(m, -np.inf)
    lw[ok] = (gammaln(m - k[ok] + 1) - gammaln(N) - gammaln(m - k[ok] - N + 2)) - (gammaln(m + 1) - gammaln(N + 1) - gammaln(m - N + 1))
    return float(np.dot(np.exp(lw), f))


def boot_med(v):
    v = np.asarray(v); b = [np.median(rng.choice(v, len(v), replace=True)) for _ in range(BOOT)]
    return float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def holm(ps):
    ps = np.asarray(ps, float); o = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0
    for i, j in enumerate(o): run = max(run, min(1.0, (m - i) * ps[j])); adj[j] = run
    return adj.tolist()


out = {"table6": {}, "table7": {}, "tests": {}, "n_inst": {}}
for K, n in SETTINGS:
    seeds = sorted(s for (k, nn, s) in rows if k == K and nn == n)
    if not seeds: continue
    out["n_inst"][f"{K}_{n}"] = len(seeds); ref = {s: REF[(K, n, s)] for s in seeds}
    t6 = {}; t7 = {}; per = {}
    for sv in SOLVERS:
        if sv not in rows[(K, n, seeds[0])]["runs"]: continue
        efforts = list(rows[(K, n, seeds[0])]["runs"][sv].keys())
        t6[sv] = {e: dict(gap=float(np.median([np.median([(x["f"] / ref[s] - 1) * 100 for x in rows[(K, n, s)]["runs"][sv][e]]) for s in seeds])),
                          t=float(np.mean([np.mean([x["t"] for x in rows[(K, n, s)]["runs"][sv][e]]) for s in seeds])),
                          runs=int(np.min([len(rows[(K, n, s)]["runs"][sv][e]) for s in seeds]))) for e in efforts}
        t7[sv] = {}; per[sv] = {}
        for T in BUD:
            best = None
            for e in efforts:
                v = []
                for s in seeds:
                    pool = rows[(K, n, s)]["runs"][sv][e]; tt = np.mean([x["t"] for x in pool]); N = int(T // tt)
                    v.append((exp_min_exact([x["f"] for x in pool], N) / ref[s] - 1) * 100)
                v = np.array(v)
                if np.isnan(v).any(): continue
                if best is None or np.median(v) < best[1]: best = (e, float(np.median(v)), v)
            if best is None: t7[sv][str(T)] = None; continue
            lo, hi = boot_med(best[2]); t7[sv][str(T)] = dict(effort=best[0], gap=best[1], lo=lo, hi=hi)
            per[sv][T] = best[2]
    out["table6"][f"{K}_{n}"] = t6; out["table7"][f"{K}_{n}"] = t7
    tests = []
    for a, b in (("SA-tuned", "bSB"), ("SA-tuned", "SA-frozen")):
        for T in BUD:
            if T in per.get(a, {}) and T in per.get(b, {}):
                dv = per[a][T] - per[b][T]
                p = float(wilcoxon(dv).pvalue) if np.any(dv != 0) else 1.0
                tests.append(dict(a=a, b=b, T=T, median_diff=float(np.median(dv)), wins=int((dv < 0).sum()), p=p))
    for fam in ({"bSB"}, {"SA-frozen"}):
        sel = [t for t in tests if t["b"] in fam]
        for t, pa in zip(sel, holm([t["p"] for t in sel])): t["p_holm"] = pa
    out["tests"][f"{K}_{n}"] = tests
json.dump(out, open("results/fresh67_summary.json", "w"), indent=1)

# Figures 1 and 2 (same layout as before: one panel per setting)
STY = {"SA-tuned": ("native SA (tuned $T_0$)", "#1b6ca8", "-"), "SA-frozen": ("native SA ($T_0 = 2$)", "#1b6ca8", "--"),
       "bSB": ("bSB (untuned)", "#c0392b", "-"), "bSB-tuned": ("bSB (tuned)", "#c0392b", "--")}
for K, sets, fn in ((3, [(3, 200), (3, 500)], "fig_K3"), (4, [(4, 200)], "fig_K4")):
    sets = [s for s in sets if f"{s[0]}_{s[1]}" in out["table7"]]
    if not sets: continue
    fig, axs = plt.subplots(1, len(sets), sharey=True, figsize=(3.45 * len(sets) if len(sets) > 1 else 4.6, 2.9 if len(sets) > 1 else 2.9), squeeze=False)
    for ax, (k, n) in zip(axs[0], sets):
        t7 = out["table7"][f"{k}_{n}"]
        for sv, d in t7.items():
            xs = [T for T in BUD if d.get(str(T))]; lab, col, ls = STY[sv]
            ax.plot(xs, [d[str(T)]["gap"] for T in xs], ls, color=col, marker="o", ms=3, lw=1.2, label=lab)
            ax.fill_between(xs, [d[str(T)]["lo"] for T in xs], [d[str(T)]["hi"] for T in xs], color=col, alpha=0.12, lw=0)
        ax.set_xscale("log"); ax.set_xlabel("Wall-clock budget (s)", fontsize=8); (ax.set_ylabel("Gap to best-known (%)", fontsize=8) if ax is axs[0][0] else None)
        ax.set_title(f"$K = {k}$, $n = {n}$", fontsize=9); ax.tick_params(labelsize=7); ax.grid(alpha=0.3, lw=0.5)
    h_, l_ = axs[0][0].get_legend_handles_labels()
    fig.legend(h_, l_, loc="upper center", ncol=4 if len(sets) > 1 else 3, fontsize=6.5, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.88 if len(sets) > 1 else 0.90)); fig.savefig(f"paper_run/figures/{fn}.png", dpi=600); fig.savefig(f"paper_run/figures/{fn}.pdf"); plt.close(fig)
print(json.dumps({k: {sv: {T: (round(v["gap"], 1) if v else None) for T, v in d.items()} for sv, d in t.items()} for k, t in out["table7"].items()}, indent=0))

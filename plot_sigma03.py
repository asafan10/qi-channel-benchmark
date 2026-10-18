"""Figure for the full refinement protocol on the irregular layout (jitter 0.3, K = 3, 200 cells, 30 instances; Table 12 and Figure 4 of the paper).
Left: exact expected best-of-N gap (median over instances, 95% bootstrap band) against wall-clock budget.
Right: median paired difference (native SA minus bSB + polish) with 95% bootstrap interval and the +-0.5-point equivalence margin.
Reads results/sigma03.json (analyse_sigma_ext.py) -> figures/sigma03_matched.png (600 dpi) and .pdf"""
import json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from figstyle import STY as _FS, XLAB, YLAB, FONT

S = json.load(open("results/sigma03.json")); B = (1, 3, 10, 30, 100)
plt.rcParams.update(FONT)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.16, 2.9))
sty = {s: _FS["SA-frozen" if s == "SA-native" else s] for s in ("SA-native", "SA-tuned", "bSB", "bSB+polish")}
for s, (c, mk, ls, lab) in sty.items():
    pts = [(T, S["res"][f"{s}|{T}"]) for T in B if S["res"].get(f"{s}|{T}")]
    x = [p[0] for p in pts]; m = [p[1][0] for p in pts]
    a1.plot(x, m, color=c, marker=mk, ls=ls, label=lab, ms=4, lw=1.2)
    a1.fill_between(x, [p[1][1] for p in pts], [p[1][2] for p in pts], color=c, alpha=.15)
a1.set_xscale("log"); a1.set_xlabel(XLAB); a1.set_ylabel(YLAB); a1.set_ylim(0, 7.6)
a1.grid(alpha=.3, which="both"); a1.legend(frameon=False, ncol=2, loc="upper center", columnspacing=0.8); a1.set_title("(a) Gap to the best plan found")
tp = {(t["a"], t["T"]): t for t in S["tests_vs_polish"]}
for k, (a, c, mk, lab, dx) in enumerate((("SA-tuned", _FS["SA-tuned"][0], _FS["SA-tuned"][1], "tuned native SA minus bSB + polish", 0.94), ("SA-native", _FS["SA-frozen"][0], _FS["SA-frozen"][1], r"native SA ($T_0 = 2$) minus bSB + polish", 1.06))):
    ts = [tp[(a, T)] for T in B if (a, T) in tp]
    x = [t["T"] * dx for t in ts]; m = [t["med_diff"] for t in ts]
    a2.errorbar(x, m, yerr=[[t["med_diff"] - t["ci"][0] for t in ts], [t["ci"][1] - t["med_diff"] for t in ts]], color=c, marker=mk, ls="none",
                ms=4, capsize=2, lw=1.0, label=lab)
a2.axhspan(-0.5, 0.5, color="grey", alpha=.12, lw=0); a2.axhline(0, color="k", lw=0.8)
for y in (-0.5, 0.5): a2.axhline(y, color="grey", lw=0.9, ls=":")
a2.text(1.0, 0.54, "equivalence margin ±0.5 points", fontsize=7.5, color="0.25")
a2.set_xscale("log"); a2.set_xlabel(XLAB); a2.set_ylabel("Median paired difference (points)"); a2.set_ylim(-0.75, 0.75)
a2.grid(alpha=.3, which="both"); a2.legend(frameon=False, loc="lower left"); a2.set_title("(b) Paired difference; negative favours native SA")
fig.tight_layout(w_pad=1.2); fig.savefig("figures/sigma03_matched.pdf"); fig.savefig("figures/sigma03_matched.png", dpi=600)
print("ok")

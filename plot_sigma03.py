"""Figure for the full refinement protocol on the irregular layout (jitter 0.3, K = 3, 200 cells, 30 instances; Table 15 of the paper).
Left: exact expected best-of-N gap (median over instances, 95% bootstrap band) against wall-clock budget.
Right: median paired difference (native SA minus bSB + polish) with 95% bootstrap interval and the +-0.5-point equivalence margin.
Reads results/sigma03.json (analyse_sigma_ext.py) -> figures/sigma03_matched.png (600 dpi) and .pdf"""
import json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = json.load(open("results/sigma03.json")); B = (1, 3, 10, 30, 100)
plt.rcParams.update({"font.size": 8})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.16, 2.9))
sty = {"SA-native": ("C1", "o", "-", "native SA (frozen $T_0$)"), "SA-tuned": ("C5", "D", "-.", "native SA (tuned $T_0$)"),
       "bSB": ("C4", "s", "--", "bSB alone"), "bSB+polish": ("C3", "^", "-", "bSB + native SA polish")}
for s, (c, mk, ls, lab) in sty.items():
    pts = [(T, S["res"][f"{s}|{T}"]) for T in B if S["res"].get(f"{s}|{T}")]
    x = [p[0] for p in pts]; m = [p[1][0] for p in pts]
    a1.plot(x, m, color=c, marker=mk, ls=ls, label=lab, ms=4, lw=1.2)
    a1.fill_between(x, [p[1][1] for p in pts], [p[1][2] for p in pts], color=c, alpha=.15)
a1.set_xscale("log"); a1.set_xlabel("wall-clock budget T (s)"); a1.set_ylabel("expected best-of-N gap (%)"); a1.set_ylim(0, 7.6)
a1.grid(alpha=.3, which="both"); a1.legend(fontsize=7, frameon=False, ncol=2, loc="upper center", columnspacing=0.8); a1.set_title("(a) gap to the best plan found", fontsize=8)
tp = {(t["a"], t["T"]): t for t in S["tests_vs_polish"]}
for k, (a, c, mk, lab, dx) in enumerate((("SA-tuned", "C5", "D", "tuned SA minus bSB + polish", 0.94), ("SA-native", "C1", "o", "frozen SA minus bSB + polish", 1.06))):
    ts = [tp[(a, T)] for T in B if (a, T) in tp]
    x = [t["T"] * dx for t in ts]; m = [t["med_diff"] for t in ts]
    a2.errorbar(x, m, yerr=[[t["med_diff"] - t["ci"][0] for t in ts], [t["ci"][1] - t["med_diff"] for t in ts]], color=c, marker=mk, ls="none",
                ms=4, capsize=2, lw=1.0, label=lab)
a2.axhspan(-0.5, 0.5, color="grey", alpha=.12, lw=0); a2.axhline(0, color="k", lw=0.8)
for y in (-0.5, 0.5): a2.axhline(y, color="grey", lw=0.9, ls=":")
a2.text(1.0, 0.54, "equivalence margin ±0.5 point", fontsize=8, color="0.25")
a2.set_xscale("log"); a2.set_xlabel("wall-clock budget T (s)"); a2.set_ylabel("median paired difference (points)"); a2.set_ylim(-0.75, 0.75)
a2.grid(alpha=.3, which="both"); a2.legend(fontsize=7, frameon=False, loc="lower left"); a2.set_title("(b) negative favours native SA", fontsize=8)
fig.tight_layout(w_pad=1.2); fig.savefig("figures/sigma03_matched.pdf"); fig.savefig("figures/sigma03_matched.png", dpi=600)
print("ok")

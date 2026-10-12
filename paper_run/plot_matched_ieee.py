"""Figures 1 and 2 of the paper (matched wall-clock, protocol P4, five instances), regenerated from the same computation as Table 7:
analyse_tuned.py (K = 3, SA-native / untuned bSB / tuned bSB) is executed unchanged for its estimates; only the figure is drawn here,
at IEEE column/page width with 8-pt text and without the in-figure footnote (the caveats are in the caption and text).
Usage (in paper_run/): python3 plot_matched_ieee.py 3|4  -> figures/fig_K3.png, figures/fig_K4.png (600 dpi) and .pdf"""
import sys, io, contextlib
src = open("analyse_tuned.py").read(); src = src[:src.index("fig, axs = plt.subplots")]
with contextlib.redirect_stdout(io.StringIO()):
    exec(src)
plt.rcParams.update({"font.size": 8})
w_ = 7.16 if len(sizes) > 1 else 3.5
fig, axs = plt.subplots(1, len(sizes), figsize=(w_, 2.7), squeeze=False)
sty = {"SA-native": ("C1", "o", "-", "native SA"), "bSB": ("C4", "s", "--", "bSB (untuned)"), "bSB-tuned": ("C2", "^", "-.", "bSB (tuned, frozen)")}
for ax, n in zip(axs[0], sizes):
    for w, _ in NAMES:
        pts = [(T, res[(n, w, T)]) for T in BUDGETS if res[(n, w, T)]]
        if not pts: continue
        x = [p[0] for p in pts]; m = [p[1][1] for p in pts]; lo = [p[1][2] for p in pts]; hi = [p[1][3] for p in pts]
        c, mk, ls, lab = sty[w]; ax.plot(x, m, color=c, marker=mk, ls=ls, label=lab, ms=4, lw=1.2); ax.fill_between(x, lo, hi, color=c, alpha=.15)
    ax.axvline(1, color="k", lw=0.8, ls="--", alpha=.5); ax.axhline(5, color="grey", lw=0.9, ls=":")
    ax.set_xscale("log"); ax.set_xlabel("wall-clock budget T (s)"); ax.set_ylim(bottom=0)
    if ax is axs[0][0]: ax.set_ylabel("expected best-of-N gap (%)")
    ax.set_title(f"K = {K}, n = {n}", fontsize=8); ax.grid(alpha=.3, which="both")
h, l = axs[0][0].get_legend_handles_labels()
if len(sizes) > 1:
    fig.legend(h, l, loc="upper center", ncol=len(h), frameon=False, bbox_to_anchor=(0.5, 1.0)); fig.tight_layout(rect=(0, 0, 1, 0.9), w_pad=1.0)
else:
    axs[0][0].legend(fontsize=7, frameon=False, loc="upper right"); fig.tight_layout()
fig.savefig(f"figures/fig_K{K}.pdf"); fig.savefig(f"figures/fig_K{K}.png", dpi=600)
for n in sizes: print(K, n, [(T, [round(res[(n, w, T)][1], 1) if res[(n, w, T)] else None for w, _ in NAMES]) for T in BUDGETS])

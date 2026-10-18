"""One visual convention for Figures 1-4: each method keeps its colour, marker, line style and name in every figure; axes and panel
titles use the symbols of the paper (wall-clock budget tau, best-of-N, K, n)."""
STY = {
    "SA-frozen": ("C1", "o", "-", r"native SA ($T_0 = 2$)"),
    "SA-tuned": ("C5", "D", "-.", r"native SA (tuned $T_0$)"),
    "bSB": ("C4", "s", "--", "bSB (untuned)"),
    "bSB-tuned": ("C4", "P", ":", "bSB (tuned)"),
    "bSB+descent": ("C0", "v", ":", "bSB + 1-opt descent"),
    "bSB+polish": ("C3", "^", "-", "bSB + native SA polish"),
}
XLAB = r"Wall-clock budget $\tau$ (s)"
YLAB = r"Expected best-of-$N$ gap (%)"
FONT = {"font.size": 8, "legend.fontsize": 7.5, "axes.titlesize": 8.5}
def title(K, n): return rf"$K = {K}$, $n = {n}$"

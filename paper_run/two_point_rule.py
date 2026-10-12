"""Basis of the 2-point reading rule for the five-instance tables: half-widths of the 95% bootstrap intervals of the
medians in Table 7 (native SA and untuned bSB, K = 3 at 200 and 500 cells and K = 4 at 200 cells, five budgets), computed by the
same code as Table 7 (analyse_tuned.py). -> ../results/two_point_rule.json"""
import io, sys, json, contextlib, numpy as np
hw = {}
for K in (3, 4):
    src = open("analyse_tuned.py").read(); src = src[:src.index("fig, axs = plt.subplots")]
    g = {"__name__": "x"}; sys.argv = ["x", str(K)]
    with contextlib.redirect_stdout(io.StringIO()): exec(src, g)
    for (n, w, T), b in g["res"].items():
        if b and w in ("SA-native", "bSB"): hw[f"{K}|{n}|{w}|{T}"] = float((b[3] - b[2]) / 2)
v = np.array(list(hw.values()))
out = dict(n_cells=len(v), median=float(np.median(v)), q25=float(np.percentile(v, 25)), q75=float(np.percentile(v, 75)), cells=hw)
json.dump(out, open("../results/two_point_rule.json", "w"), indent=1); print({k: out[k] for k in ("n_cells", "median", "q25", "q75")})

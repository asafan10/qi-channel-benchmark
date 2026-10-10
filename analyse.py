import json, sys, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

tag = sys.argv[1]; files = sys.argv[2:]
merged = {}
for f in files:                       # later files override solver entries of earlier ones
    for l in open(f):
        r = json.loads(l); key = (r["n"], r["seed"])
        if key in merged: merged[key]["res"].update(r["res"])
        else: merged[key] = r
rows = list(merged.values())
order = ["Greedy", "Tabu", "SA-native", "CP-SAT", "SA-QUBO", "PT-QUBO", "bSB", "dSB"]
sizes = sorted({r["n"] for r in rows})
out = Path("figures"); out.mkdir(exist_ok=True)

gap = {s: {n: [] for n in sizes} for s in order}
tim = {s: {n: [] for n in sizes} for s in order}
wins = {s: {n: 0 for n in sizes} for s in order}
opt_cnt = {n: [0, 0] for n in sizes}
for r in rows:
    n = r["n"]; res = r["res"]
    bks = min(v["E"] for v in res.values())
    for s in order:
        g = 100 * (res[s]["E"] - bks) / max(bks, 1e-9)
        gap[s][n].append(g); tim[s][n].append(res[s]["t"])
        if res[s]["E"] <= bks * (1 + 1e-6): wins[s][n] += 1
    opt_cnt[n][0] += int(res["CP-SAT"].get("optimal", False)); opt_cnt[n][1] += 1

lines = ["| n cells (vars) | " + " | ".join(order) + " |", "|---|" + "---|" * len(order)]
for n in sizes:
    cells = [f"{n} ({3*n})"]
    for s in order:
        cells.append(f"{np.mean(gap[s][n]):.1f} ± {np.std(gap[s][n]):.1f} [{wins[s][n]}/{len(gap[s][n])}]")
    lines.append("| " + " | ".join(cells) + " |")
tab = "\n".join(lines)
lines2 = ["| n | " + " | ".join(order) + " |", "|---|" + "---|" * len(order)]
for n in sizes:
    lines2.append(f"| {n} | " + " | ".join(f"{np.mean(tim[s][n]):.2f}" for s in order) + " |")
tab2 = "\n".join(lines2)
Path(f"table_{tag}.md").write_text(
    "Gap to best-known (%), mean ± std over instances [times matching best-known]\n\n" + tab +
    "\n\nMean wall-clock per solver run (s)\n\n" + tab2 +
    "\n\nCP-SAT proven optimal: " + ", ".join(f"n={n}: {opt_cnt[n][0]}/{opt_cnt[n][1]}" for n in sizes) + "\n")
print(Path(f"table_{tag}.md").read_text())

# Fig A: gap vs size
fig, ax = plt.subplots(figsize=(6.4, 4))
mk = dict(zip(order, "os^Dvp*X"))
for s in order:
    ax.plot(sizes, [np.mean(gap[s][n]) for n in sizes], marker=mk[s], label=s)
ax.set_xscale("log"); ax.set_xlabel("cells n (binary vars = 3n)"); ax.set_ylabel("gap to best-known (%)")
ax.set_title(f"Solution quality vs scale ({tag}); mean of instances"); ax.grid(alpha=.3); ax.legend(ncol=2, fontsize=8)
fig.tight_layout(); fig.savefig(out / f"gap_vs_size_{tag}.pdf")

# Fig B: convergence at largest n
n = sizes[-1]
fig, ax = plt.subplots(figsize=(6.4, 4))
sel = [r for r in rows if r["n"] == n]
tmax = max(v["t"] for r in sel for v in r["res"].values()) 
grid = np.geomspace(0.05, tmax, 60)
for s in order:
    if s == "Greedy": continue
    curves = []
    for r in sel:
        c = r["res"][s]["curve"]; ts = np.array([p[0] for p in c]); es = np.array([p[1] for p in c])
        bks = min(v["E"] for v in r["res"].values())
        y = np.array([es[ts <= t].min() if (ts <= t).any() else np.nan for t in grid])
        curves.append(100 * (y - bks) / bks)
    ax.plot(grid, np.nanmedian(np.array(curves), axis=0), label=s)
ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=1); ax.set_xlabel("wall-clock (s)")
ax.set_ylabel("median gap to best-known (%)"); ax.set_title(f"Anytime convergence, n={n} cells ({tag})")
ax.grid(alpha=.3); ax.legend(ncol=2, fontsize=8); fig.tight_layout(); fig.savefig(out / f"convergence_{tag}.pdf")
print("figures written")

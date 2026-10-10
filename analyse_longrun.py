import json, sys, numpy as np
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

K = int(sys.argv[1]) if len(sys.argv) > 1 else 3
BUDGETS = [1, 3, 10, 30, 100]
rng = np.random.default_rng(7)


def jl(f): return [json.loads(l) for l in open(f)] if Path(f).exists() else []


sweep = {(r["n"], r["seed"]): r for r in jl(f"sweep_K{K}.jsonl")}
long_ = {(r["n"], r["seed"]): r for r in jl(f"longrun_K{K}.jsonl")}
ref = {}
for f in (f"sweep_K{K}.jsonl", f"tts_K{K}.jsonl", f"refrun_K{K}.jsonl", f"longrun_K{K}.jsonl", f"tuned_K{K}.jsonl"):
    for r in jl(f):
        k = (r["n"], r["seed"]); ref[k] = min(ref.get(k, np.inf), r.get("prev_bks", np.inf), r["new_best"])

sizes = sorted({k[0] for k in long_}); out = []


def samples(n, seed, solver):
    """-> dict effort factor -> (np.array of f per run, mean time per run, n_runs); bSB also replica-level."""
    d = {}; rep = {}
    if (n, seed) in sweep:
        for fk, c in sweep[(n, seed)]["runs"][solver].items():
            if solver in ("SA-native", "bSB"): d[float(fk)] = (np.array(c["f"]), float(np.mean(c["t"])))
    if (n, seed) in long_:
        for fk, recs in long_[(n, seed)]["runs"][solver].items():
            d[float(fk)] = (np.array([r["f_sel"] for r in recs]), float(np.mean([r["t"] for r in recs])))
            if solver == "bSB":
                rep[float(fk)] = (np.concatenate([r["f_all"] for r in recs]), float(np.mean([r["t"] for r in recs])) / 8.0)
    return d, rep


def exp_min(f, N, B=400):
    if N >= 1000: N = 1000
    return float(np.mean([rng.choice(f, size=N, replace=True).min() for _ in range(B)]))


def best_at_budget(n, seeds, solver, T, replica=False):
    facts = None; per = {}
    for s in seeds:
        d, rep = samples(n, s, solver); src = rep if replica else d
        for fk, (f, t) in src.items():
            N = int(T // t)
            if N >= 1 and len(f) >= 1:
                per.setdefault(fk, []).append((exp_min(f, N) / ref[(n, s)] - 1) * 100)
            else:
                per.setdefault(fk, []).append(np.nan)
    best = None
    for fk, v in per.items():
        if len(v) == len(seeds) and not np.any(np.isnan(v)):
            m = np.median(v)
            if best is None or m < best[1]: best = (fk, m)
    return best  # (factor, median gap %)


lines = [f"# Symmetric long-run comparison, K={K}\n"]
curves = {}
for n in sizes:
    seeds = sorted(s for (nn, s) in long_ if nn == n)
    lines.append(f"\n## n = {n} cells ({len(seeds)} instances)\n")
    lines.append("### Single-run median gap to updated best-known (%) and mean time per run (ms)\n")
    lines.append("| solver | effort factor | runs/instance | median gap % | time/run (ms) |"); lines.append("|---|---|---|---|---|")
    for solver in ("SA-native", "bSB"):
        allf = sorted({fk for s in seeds for fk in samples(n, s, solver)[0]})
        for fk in allf:
            gaps = []; ts = []; nr = []
            for s in seeds:
                d, _ = samples(n, s, solver)
                if fk in d: gaps.append(np.median((d[fk][0] / ref[(n, s)] - 1) * 100)); ts.append(d[fk][1]); nr.append(len(d[fk][0]))
            if len(gaps) == len(seeds):
                lines.append(f"| {solver} | x{fk:.3g} | {int(np.median(nr))} | {np.median(gaps):.1f} | {1000*np.mean(ts):.0f} |")
    lines.append("\n### Matched wall-clock: expected best-of-N gap (%) with the best effort level chosen at each budget (median over instances)\n")
    lines.append("| budget T | SA-native | bSB (batch = 1 run) | bSB (each replica = 1 run, optimistic) |"); lines.append("|---|---|---|---|")
    curves[n] = {"SA-native": [], "bSB": [], "bSB-rep": []}
    for T in BUDGETS:
        row = []
        for key, solver, rep in (("SA-native", "SA-native", False), ("bSB", "bSB", False), ("bSB-rep", "bSB", True)):
            b = best_at_budget(n, seeds, solver, T, rep)
            curves[n][key].append(b[1] if b else np.nan)
            row.append(f"{b[1]:.1f}% (x{b[0]:.3g})" if b else "-")
        lines.append(f"| {T} s | " + " | ".join(row) + " |")
txt = "\n".join(lines) + "\n"; Path(f"longrun_K{K}.md").write_text(txt); print(txt)

fig, axs = plt.subplots(1, len(sizes), figsize=(5.2 * len(sizes), 4), squeeze=False)
for ax, n in zip(axs[0], sizes):
    for key, lab, c in (("SA-native", "SA-native", "C1"), ("bSB", "bSB (batch = run)", "C4"), ("bSB-rep", "bSB (replica = run)", "C4")):
        ax.plot(BUDGETS, curves[n][key], marker="o", label=lab, color=c, ls="--" if key == "bSB-rep" else "-")
    ax.set_xscale("log"); ax.set_xlabel("wall-clock budget T (s)"); ax.set_ylabel("expected best-of-N gap to best-known (%)")
    ax.set_title(f"K={K}, n={n}: matched wall-clock", fontsize=9); ax.grid(alpha=.3, which="both"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"figures/longrun_K{K}.pdf")

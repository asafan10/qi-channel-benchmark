import json, sys, numpy as np
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

K = int(sys.argv[1]) if len(sys.argv) > 1 else 3; BUDGETS = [1, 3, 10, 30, 100]; rng = np.random.default_rng(11); BOOT = 2000


def jl(f): return [json.loads(l) for l in open(f)] if Path(f).exists() else []


sweep = {(r["n"], r["seed"]): r for r in jl(f"sweep_K{K}.jsonl")}
long_ = {(r["n"], r["seed"]): r for r in jl(f"longrun_K{K}.jsonl")}
tuned = {(r["n"], r["seed"]): r for r in jl(f"tuned_K{K}.jsonl")}
ref = {}
for f in (f"sweep_K{K}.jsonl", f"tts_K{K}.jsonl", f"refrun_K{K}.jsonl", f"longrun_K{K}.jsonl", f"tuned_K{K}.jsonl"):
    for r in jl(f):
        k = (r["n"], r["seed"]); ref[k] = min(ref.get(k, np.inf), r.get("prev_bks", np.inf), r["new_best"])


def samples(n, s, which):
    d = {}
    if which in ("SA-native", "bSB"):
        if (n, s) in sweep:
            for fk, c in sweep[(n, s)]["runs"][which].items(): d[float(fk)] = (np.array(c["f"]), float(np.mean(c["t"])))
        if (n, s) in long_:
            for fk, recs in long_[(n, s)]["runs"][which].items():
                d[float(fk)] = (np.array([r["f_sel"] for r in recs]), float(np.mean([r["t"] for r in recs])))
    else:  # bSB-tuned
        if (n, s) in tuned:
            for fk, recs in tuned[(n, s)]["runs"].items():
                d[float(fk)] = (np.array([r["f_sel"] for r in recs]), float(np.mean([r["t"] for r in recs])))
    return d


def exp_min(f, N, B=300):
    N = min(N, 1000); return float(np.mean([rng.choice(f, size=N, replace=True).min() for _ in range(B)]))


def at_budget(n, seeds, which, T):
    per = {}
    for s in seeds:
        for fk, (f, t) in samples(n, s, which).items():
            N = int(T // t)
            per.setdefault(fk, {})[s] = (exp_min(f, N) / ref[(n, s)] - 1) * 100 if N >= 1 else np.nan
    best = None
    for fk, dct in per.items():
        v = np.array([dct.get(s, np.nan) for s in seeds])
        if np.isnan(v).any(): continue
        m = np.median(v)
        if best is None or m < best[1]: best = (fk, m, v)
    if best is None: return None
    v = best[2]; bs = [np.median(rng.choice(v, size=len(v), replace=True)) for _ in range(BOOT)]
    return best[0], best[1], np.percentile(bs, 2.5), np.percentile(bs, 97.5)


NAMES = [("SA-native", "SA-native"), ("bSB", "bSB (untuned)"), ("bSB-tuned", "bSB (tuned, frozen)")]
BASE = tuned if tuned else long_
NAMES = [x for x in NAMES if x[0] != "bSB-tuned" or tuned]
sizes = sorted({k[0] for k in BASE}); lines = [f"# Matched wall-clock, K={K}: SA-native vs untuned vs tuned bSB (median over instances; [95% bootstrap CI over instances]; best effort factor)\n"]
res = {}
for n in sizes:
    seeds = sorted(s for (nn, s) in BASE if nn == n)
    lines += [f"\n## n = {n} ({len(seeds)} instances)\n", "| T | " + " | ".join(l for _, l in NAMES) + " |", "|---|" + "---|" * len(NAMES)]
    for T in BUDGETS:
        row = []
        for w, _ in NAMES:
            b = at_budget(n, seeds, w, T); res[(n, w, T)] = b
            row.append(f"{b[1]:.1f}% [{b[2]:.1f}, {b[3]:.1f}] (x{b[0]:.3g})" if b else "-")
        lines.append(f"| {T} s | " + " | ".join(row) + " |")
    lines += [f"\n### single-run median gap (%) / time per run (ms), tuned bSB, n={n}\n"] if tuned else []
    for fk in (sorted({fk for s in seeds for fk in samples(n, s, "bSB-tuned")}) if tuned else []):
        g = np.median([np.median((samples(n, s, "bSB-tuned")[fk][0] / ref[(n, s)] - 1) * 100) for s in seeds])
        t = np.mean([samples(n, s, "bSB-tuned")[fk][1] for s in seeds]); lines.append(f"- x{fk:.3g}: {g:.1f}% ({1000*t:.0f} ms)")
txt = "\n".join(lines) + "\n"; Path(f"tuned_K{K}.md").write_text(txt); print(txt)

fig, axs = plt.subplots(1, len(sizes), figsize=(5.4 * len(sizes), 4.2), squeeze=False)
sty = {"SA-native": ("C1", "o", "-"), "bSB": ("C4", "s", "--"), "bSB-tuned": ("C2", "^", "-")}
for ax, n in zip(axs[0], sizes):
    for w, lab in NAMES:
        pts = [(T, res[(n, w, T)]) for T in BUDGETS if res[(n, w, T)]]
        if not pts: continue
        x = [p[0] for p in pts]; m = [p[1][1] for p in pts]; lo = [p[1][2] for p in pts]; hi = [p[1][3] for p in pts]
        c, mk, ls = sty[w]; ax.plot(x, m, color=c, marker=mk, ls=ls, label=lab); ax.fill_between(x, lo, hi, color=c, alpha=.15)
    ax.axvline(1, color="k", lw=1.0, ls="--", alpha=.5)
    ax.axhline(5, color="grey", lw=1.2, ls=":"); ax.text(ax.get_xlim()[1] * 0.97, 5.2, "5% gap", fontsize=9, ha="right", color="grey")
    ax.set_xscale("log"); ax.set_xlabel("wall-clock budget T (s)"); ax.set_ylabel("expected best-of-N gap to best-known (%)")
    ax.set_title(f"K={K}, n={n}", fontsize=10); ax.grid(alpha=.3, which="both"); ax.legend(fontsize=8)
fig.text(0.5, -0.07, "Best-known is heuristic and includes the SA-native runs themselves, so SA-native values near 0% at large T are reference-limited.\nBands: 95% bootstrap CI over 5 instances. Dotted line: 5% gap. Dashed vertical line: 1 s budget.", ha="center", fontsize=10)
fig.tight_layout(); fig.savefig(f"figures/final_matched_wallclock_K{K}.pdf", bbox_inches="tight")

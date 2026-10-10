import json, sys, numpy as np
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

SOLVERS = ["Tabu", "SA-native", "SA-QUBO", "PT-QUBO", "bSB", "dSB"]
COL = {"Tabu": "C0", "SA-native": "C1", "SA-QUBO": "C2", "PT-QUBO": "C3", "bSB": "C4", "dSB": "C5"}
EPS = [1.0, 2.0, 5.0]
B = 2000
rng = np.random.default_rng(12345)


def tts_from(p, trun, conf=0.99):
    if p <= 0: return np.inf
    if p >= conf: return trun
    return trun * np.log(1 - conf) / np.log(1 - p)


def load(f):
    rows = [json.loads(l) for l in open(f)]
    for r in rows:
        r["ref"] = min(r["prev_bks"], r["new_best"])
    return rows


def inst_stats(r, s, eps):
    f = np.array(r["runs"][s]["f"]); t = np.array(r["runs"][s]["t"])
    return (f <= r["ref"] * (1 + eps / 100) + 1e-9), t


def analyse(rows, K):
    sizes = sorted({r["n"] for r in rows}); res = {}
    lines = []
    for eps in EPS:
        lines.append(f"\n### K={K}, epsilon = {eps:.0f}%: median TTS(99%) in ms over instances [95% bootstrap CI]; pooled per-run success p\n")
        lines.append("| n | " + " | ".join(SOLVERS) + " |"); lines.append("|---|" + "---|" * len(SOLVERS))
        for n in sizes:
            sel = [r for r in rows if r["n"] == n]; cells = []
            for s in SOLVERS:
                ok = []; tr = []
                for r in sel:
                    o, t = inst_stats(r, s, eps); ok.append(o); tr.append(t)
                tts0 = [tts_from(o.mean(), t.mean()) for o, t in zip(ok, tr)]
                med = np.median(tts0)
                bs = np.empty(B)
                for b in range(B):
                    idx = rng.integers(0, len(sel), len(sel)); v = []
                    for i in idx:
                        rr = rng.integers(0, len(ok[i]), len(ok[i]))
                        v.append(tts_from(ok[i][rr].mean(), tr[i][rr].mean()))
                    bs[b] = np.median(v)
                lo, hi = np.percentile(bs, [2.5, 97.5])
                pooled = np.concatenate(ok).mean()
                res[(eps, n, s)] = (med, lo, hi, pooled)
                f = lambda x: "inf" if not np.isfinite(x) else f"{1000*x:.0f}"
                cells.append(f"{f(med)} [{f(lo)}, {f(hi)}]; p={pooled:.3f}")
            lines.append(f"| {n} | " + " | ".join(cells) + " |")
    txt = "\n".join(lines) + "\n"
    Path(f"tts_K{K}.md").write_text(txt); print(txt)
    return sizes, res


def fig(sizes, res, K, eps=2.0):
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    for s in SOLVERS:
        xs, m, lo, hi = [], [], [], []
        for n in sizes:
            med, l, h, p = res[(eps, n, s)]
            if np.isfinite(med): xs.append(n); m.append(1000 * med); lo.append(1000 * l); hi.append(1000 * h if np.isfinite(h) else np.nan)
        if not xs: continue
        m, lo, hi = map(np.array, (m, lo, hi))
        err = np.vstack([m - lo, np.where(np.isfinite(hi), hi - m, 0)])
        ax.errorbar(xs, m, yerr=err, marker="o", capsize=3, label=s, color=COL[s])
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("cells n"); ax.set_ylabel("TTS$_{99}$ to within %.0f%% of best-known (ms)" % eps)
    ax.set_title(f"Independent-run TTS, K={K} (median of 5 instances, 95% bootstrap CI)", fontsize=9)
    ax.grid(alpha=.3, which="both"); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(f"figures/tts_K{K}_eps{int(eps)}.pdf")


if __name__ == "__main__":
    for K in (3, 4):
        f = f"tts_K{K}.jsonl"
        if not Path(f).exists(): continue
        rows = load(f); sizes, res = analyse(rows, K)
        for e in (2.0, 5.0): fig(sizes, res, K, e)

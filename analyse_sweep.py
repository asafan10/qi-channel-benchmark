import json, sys, numpy as np
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

SOLVERS = ["Tabu", "SA-native", "SA-QUBO", "PT-QUBO", "bSB", "dSB"]
COL = {"Tabu": "C0", "SA-native": "C1", "SA-QUBO": "C2", "PT-QUBO": "C3", "bSB": "C4", "dSB": "C5"}
EPS = [5.0, 2.0]
B = 2000
rng = np.random.default_rng(2024)


def tts_from(p, trun, conf=0.99):
    if p <= 0: return np.inf
    if p >= conf: return trun
    return trun * np.log(1 - conf) / np.log(1 - p)


def load_refs(K):
    """best-known reference = min over every source we have for that instance."""
    ref = {}
    for f in (f"sweep_K{K}.jsonl", f"tts_K{K}.jsonl", f"refrun_K{K}.jsonl"):
        if not Path(f).exists(): continue
        for l in open(f):
            r = json.loads(l); k = (r["n"], r["seed"])
            ref[k] = min(ref.get(k, np.inf), r["prev_bks"], r["new_best"])
    return ref


def analyse(K):
    f = f"sweep_K{K}.jsonl"
    if not Path(f).exists(): return
    rows = [json.loads(l) for l in open(f)]; ref = load_refs(K)
    sizes = sorted({r["n"] for r in rows}); lines = []; store = {}
    for eps in EPS:
        lines.append(f"\n### K={K}, eps={eps:.0f}%: run-length-optimised TTS99 in ms, median over instances [95% bootstrap CI incl. re-selection of run length]; "
                     "(best effort factor, per-run success p, per-run time ms)\n")
        lines.append("| n | " + " | ".join(SOLVERS) + " |"); lines.append("|---|" + "---|" * len(SOLVERS))
        for n in sizes:
            sel = [r for r in rows if r["n"] == n]; cells = []
            for s in SOLVERS:
                facts = list(sel[0]["runs"][s].keys())
                data = {}
                for fk in facts:
                    data[fk] = [(np.array(r["runs"][s][fk]["f"]) <= ref[(r["n"], r["seed"])] * (1 + eps / 100) + 1e-9,
                                 np.array(r["runs"][s][fk]["t"])) for r in sel]
                med = {fk: np.median([tts_from(o.mean(), t.mean()) for o, t in data[fk]]) for fk in facts}
                best = min(facts, key=lambda k: med[k])
                # bootstrap with re-selection of the run length
                bs = np.empty(B)
                for b in range(B):
                    idx = rng.integers(0, len(sel), len(sel)); m = {}
                    for fk in facts:
                        v = []
                        for i in idx:
                            o, t = data[fk][i]; rr = rng.integers(0, len(o), len(o))
                            v.append(tts_from(o[rr].mean(), t[rr].mean()))
                        m[fk] = np.median(v)
                    bs[b] = min(m.values())
                lo, hi = np.percentile(bs, [2.5, 97.5])
                o_all = np.concatenate([o for o, t in data[best]]); t_all = np.concatenate([t for o, t in data[best]])
                store[(eps, n, s)] = dict(med=med, best=best, lo=lo, hi=hi, facts=facts, data=data)
                g = lambda x: "inf" if not np.isfinite(x) else f"{1000*x:.0f}"
                cells.append(f"{g(med[best])} [{g(lo)}, {g(hi)}] (x{float(best):.2g}, p={o_all.mean():.3f}, {1000*t_all.mean():.0f} ms)")
            lines.append(f"| {n} | " + " | ".join(cells) + " |")
    # single-run quality vs effort
    lines.append(f"\n### K={K}: median single-run gap to best-known (%) by effort factor (per-run time in ms)\n")
    for n in sizes:
        sel = [r for r in rows if r["n"] == n]; lines.append(f"\n**n={n}**\n"); lines.append("| solver | " + " | ".join(f"x{k}" for k in sel[0]["runs"]["Tabu"].keys()) + " |")
        lines.append("|---|---|---|---|---|")
        for s in SOLVERS:
            cells = []
            for fk in sel[0]["runs"][s].keys():
                g = np.median([np.median(100 * (np.array(r["runs"][s][fk]["f"]) - ref[(r["n"], r["seed"])]) / ref[(r["n"], r["seed"])]) for r in sel])
                t = np.mean([np.mean(r["runs"][s][fk]["t"]) for r in sel])
                cells.append(f"{g:.1f}% ({1000*t:.0f} ms)")
            lines.append(f"| {s} (factors {','.join(f'{float(k):.2g}' for k in sel[0]['runs'][s].keys())}) | " + " | ".join(cells) + " |")
    txt = "\n".join(lines) + "\n"; Path(f"sweep_K{K}.md").write_text(txt); print(txt)
    # figure: TTS vs per-run time at the largest size
    for eps in EPS:
        n = sizes[-1]; fig, ax = plt.subplots(figsize=(6.6, 4.2))
        for s in SOLVERS:
            st = store[(eps, n, s)]; xs, ys = [], []
            for fk in st["facts"]:
                t = np.mean(np.concatenate([tt for o, tt in st["data"][fk]]))
                if np.isfinite(st["med"][fk]): xs.append(1000 * t); ys.append(1000 * st["med"][fk])
            if xs: ax.plot(xs, ys, marker="o", label=s, color=COL[s])
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("per-run time (ms)"); ax.set_ylabel(f"TTS$_{{99}}$ within {eps:.0f}% (ms)")
        ax.set_title(f"Run-length sweep, K={K}, n={n} (median of 5 instances; missing = no success)", fontsize=9)
        ax.grid(alpha=.3, which="both"); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(f"figures/sweep_K{K}_n{n}_eps{int(eps)}.pdf"); plt.close(fig)


if __name__ == "__main__":
    for K in (3, 4): analyse(K)

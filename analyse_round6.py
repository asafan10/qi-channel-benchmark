"""Summary of bench_round6.py (K2, DW) and the phase-domain analysis of the K = 3 plans (PH).
Phase-domain analysis: on the unjittered lattice the proper 3-colouring is unique up to a permutation of the channels, so every
cell whose closed first-ring neighbourhood matches one permutation of the reuse pattern exactly is a 'bulk' cell of that phase;
the rest are 'wall' cells. Domains are connected groups of bulk cells (lattice adjacency) in the same phase.
-> results/round6_summary.json"""
import json, math, itertools, numpy as np
from bench_ran_qubo import interference, csr_from_W, make_instance
from bench_jitter import make_layout
from bench_reuse_polish import best_reuse

n = 200; PERMS = list(itertools.permutations(range(3)))
cols = int(math.ceil(math.sqrt(n * 1.15))); pts = []; r = 0
while len(pts) < n:
    for c in range(cols): pts.append((c + 0.5 * (r % 2), r * math.sqrt(3) / 2))
    r += 1
base = np.array(pts[:n]); D = np.linalg.norm(base[:, None] - base[None], axis=2)
NB = [np.where((D[i] > 0) & (D[i] < 1.01))[0] for i in range(n)]


def phases(c, P):
    ph = np.full(n, -1)
    for i in range(n):
        idx = np.r_[i, NB[i]]
        for k, p in enumerate(PERMS):
            if all(c[j] == p[P[j]] for j in idx): ph[i] = k; break
    return ph


def domains(ph, min_size=5):
    seen = np.zeros(n, bool); sizes = []
    for i in range(n):
        if ph[i] < 0 or seen[i]: continue
        st = [i]; seen[i] = True; s = 0
        while st:
            u = st.pop(); s += 1
            for v in NB[u]:
                if not seen[v] and ph[v] == ph[i]: seen[v] = True; st.append(v)
        sizes.append(s)
    return sum(1 for s in sizes if s >= min_size), (max(sizes) if sizes else 0)


def wall_cost(c, W, ph):
    ii, jj = np.nonzero(np.triu(W, 1)); tot = 0.0; wall = 0.0
    for i, j in zip(ii, jj):
        if c[i] == c[j]:
            tot += W[i, j]
            if ph[i] < 0 or ph[j] < 0: wall += W[i, j]
    return (wall / tot if tot > 0 else 0.0), wall


out = {}
K2 = json.load(open("results/round6_K2.json"))["rows"]; out["K2"] = {}
for sig in (0.15, 0.3):
    rr = [r for r in K2 if r["sigma"] == sig]
    g = lambda r, d, f: (np.median(r[d][f]["f"]) / r["best"] - 1) * 100
    out["K2"][str(sig)] = dict(
        bsb={f: dict(median_gap=float(np.median([g(r, "bsb", f) for r in rr])), t=float(np.mean([r["bsb"][f]["t"] for r in rr]))) for f in ("1", "4", "16")},
        sa={f: dict(median_gap=float(np.median([g(r, "sa", f) for r in rr])), t=float(np.mean([r["sa"][f]["t"] for r in rr]))) for f in ("4", "16", "64", "256")})
K2B = {(r["sigma"], r["seed"]): r["lb"] for r in json.load(open("results/round6b_K2B.json"))["rows"]}
for sig in (0.15, 0.3):
    rr = [r for r in K2 if r["sigma"] == sig]
    out["K2"][str(sig)]["best_to_bound"] = float(np.median([(r["best"] / K2B[(sig, r["seed"])] - 1) * 100 for r in rr]))
    out["K2"][str(sig)]["bsb_to_bound"] = {f: float(np.median([(np.median(r["bsb"][f]["f"]) / K2B[(sig, r["seed"])] - 1) * 100 for r in rr])) for f in ("1", "4", "16")}
    out["K2"][str(sig)]["bsb_excess_per_cell"] = {f: float(np.median([(np.median(r["bsb"][f]["f"]) - r["best"]) / 200 for r in rr])) for f in ("1", "4", "16")}
    out["K2"][str(sig)]["best_per_cell"] = float(np.median([r["best"] / 200 for r in rr]))
DW = json.load(open("results/round6_DW.json")); out["DW"] = {"rho_dw": DW["tune"]}
for K, sig in ((3, 0.15), (3, 0.3), (4, 0.15)):
    rr = [r for r in DW["rows"] if r["K"] == K and r["sigma"] == sig]
    g = lambda r, d, f: (np.median(r[d][f]["f"]) / r["ref"] - 1) * 100
    out["DW"][f"K{K}_s{sig}"] = dict(
        dw={f: dict(median_gap=float(np.median([g(r, "dw", f) for r in rr])), valid=float(np.mean([np.mean(r["dw"][f]["valid"]) for r in rr])),
                    t=float(np.mean([r["dw"][f]["t"] for r in rr]))) for f in ("1", "4", "16")},
        oh={f: dict(median_gap=float(np.median([g(r, "oh", f) for r in rr])), t=float(np.mean([r["oh"][f]["t"] for r in rr]))) for f in ("1", "4", "16")},
        sa={f: dict(median_gap=float(np.median([g(r, "sa", f) for r in rr])), t=float(np.mean([r["sa"][f]["t"] for r in rr]))) for f in ("4", "16", "64", "256")})
    out["DW"][f"K{K}_s{sig}"]["oh_excess_per_cell"] = {f: float(np.median([(np.median(r["oh"][f]["f"]) - r["ref"]) / 200 for r in rr])) for f in ("1", "4", "16")}
    out["DW"][f"K{K}_s{sig}"]["ref_per_cell"] = float(np.median([r["ref"] / 200 for r in rr]))
MG = {r["seed"]: r["plans"] for r in json.load(open("results/round6b_MG.json"))["rows"]}
out["PH"] = {}
for sig in (0.15, 0.3):
    rr = [r for r in DW["rows"] if r["K"] == 3 and r["sigma"] == sig]; res = {}
    for key in (("oh1", "oh16", "dw16", "sa1", "sa2", "sa4", "sa64", "sa256", "reuse") if sig == 0.15 else ("oh1", "oh16", "dw16", "sa64", "sa256", "reuse")):
        m = dict(gap=[], bulk=[], domains=[], largest=[], wall_share=[], wall_abs=[], hamming=[])
        for r in rr:
            W = make_instance(n, r["seed"], 3) if sig == 0.15 else make_layout(n, r["seed"], 3, sig); csr = csr_from_W(W)
            P = best_reuse(n, 3, csr)
            plans = [P] if key == "reuse" else [np.array(p) for p in (MG[r["seed"]][key] if key in ("sa1", "sa2", "sa4") else r["plans"][key])]
            for c in plans:
                ph = phases(c, P); nd, lg = domains(ph)
                m["gap"].append((interference(np.asarray(c, np.int64), *csr) / r["ref"] - 1) * 100); m["bulk"].append(float((ph >= 0).mean()))
                m["domains"].append(nd); m["largest"].append(lg); _ws, _wa = wall_cost(c, W, ph); m["wall_share"].append(_ws); m["wall_abs"].append(_wa)
                m["hamming"].append(float(min(np.mean(c != np.array(p)[P]) for p in PERMS)))
        res[key] = {k: float(np.median(v)) for k, v in m.items()}
    out["PH"][str(sig)] = res
json.dump(out, open("results/round6_summary.json", "w"), indent=1)
print(json.dumps(out, indent=1))

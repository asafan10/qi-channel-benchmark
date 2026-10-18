"""Table 8 (protocol P6) scored against the reference built without native SA's runs and without the
polish (refs.build_ref(no_native_sa=True)), same pools, efforts and estimator as analyse_hybrid.py -> results/table8_nosa.json"""
import io, json, contextlib, numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    exec(open("analyse_hybrid.py").read())
out = {}
for (K, n, a_, T), v in INST.items():
    seeds = v[1]; g = np.array(v[2])
    c = np.array([ref[(K, n, s)] * (1 + x / 100) for s, x in zip(seeds, g)]); rn = np.array([ref_nosa[(K, n, s)] for s in seeds])
    out[f"{K}|{n}|{a_}|{T}"] = dict(main=float(np.median(g)), nosa=float(np.median((c / rn - 1) * 100)), n=len(seeds))
json.dump(out, open("results/table8_nosa.json", "w"), indent=1)
for k in sorted(out): print(k, round(out[k]["main"], 2), round(out[k]["nosa"], 2), out[k]["n"])

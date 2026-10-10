# Stronger best-known references: very long native-space anneals (64x baseline sweeps), 5 runs per instance.
import json, sys, numpy as np, os
from bench_ran_qubo import *
from bench_tts import warmup, previous_bks
K=int(sys.argv[1]); sizes=[int(x) for x in sys.argv[2].split(",")]; outf=sys.argv[3]
warmup(K,CFG); prev=previous_bks(K)
with open(outf,"a") as fh:
    for n in sizes:
        for seed in range(5):
            W=make_instance(n,seed,K); csr=csr_from_W(W); rng=np.random.default_rng(777+seed)
            best=min(sa_native(rng.integers(0,K,n),K,*csr,2.0,0.02,64*CFG["native_sweeps"],900000+k)[1] for k in range(5))
            fh.write(json.dumps(dict(n=n,seed=seed,K=K,prev_bks=prev[(n,seed)],new_best=best))+"\n"); fh.flush()
            print(n,seed,"prev %.2f longSA %.2f"%(prev[(n,seed)],best),flush=True)

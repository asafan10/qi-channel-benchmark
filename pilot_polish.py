# Pilot on HELD-OUT seeds 100,101 (not evaluation seeds 0-4): choose polish temperature
import time, numpy as np
from exp_common import *
warm(3); K=3
for n in (200,500):
  for seed in (100,101):
    W, csr, a, J, lam, jj = setup(n, seed, K); rng=np.random.default_rng(5)
    Xb=sb_batch(a,J,4*3000,8,rng,False); es=[qubo_energy(x,a,*jj) for x in Xb]; rr=repair_all(Xb,K,csr); ch0=rr[int(np.argmin(es))][0]
    print(n,seed,"bSB",round(rr[int(np.argmin(es))][1],2))
    for sw in (1000,4000):
      for T0 in (2.0,0.6,0.3,0.15,0.08):
        vb=[sa_native(ch0,K,*csr,T0,0.02,sw,s)[1] for s in range(3)]
        vr=[sa_native(rng.integers(0,K,n),K,*csr,T0,0.02,sw,s)[1] for s in range(3)]
        print(f"  sweeps {sw} T0 {T0}: from bSB {np.mean(vb):.2f}  from random {np.mean(vr):.2f}")

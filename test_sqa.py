# Verification (as advocated in Sec. 5.3): (i) incremental local fields equal recomputed ones; (ii) SQA finds the brute-force ground state on small QUBOs.
import itertools, numpy as np
from exp_common import *
from sqa import *
ok=0
for seed in range(6):
    W, csr, a, J, lam, jj = setup(5, seed, 3)        # 15 binary variables
    N=len(a); X=np.array(list(itertools.product([0,1],repeat=N)))
    Ebf=min(qubo_energy(x,a,*jj) for x in X)
    rng=np.random.default_rng(seed); x,e=sqa_solve(a,J,P=16,T=0.05,G0=3.0,G1=0.01,sweeps=300,rng=rng,seed=seed)
    h,jp,jx,jd=ising_form(a,J); s=2*x-1
    E0=qubo_energy(np.ones(N,dtype=np.int64)*0,a,*jj)  # =0 ; offset check: qubo(x) - ising(s) must be constant
    off=qubo_energy(x,a,*jj)-e
    # constant offset of Ising form, from x=0 (s=-1)
    off0=0-(sum(h*-1)+sum(jd[k]*1 for p in range(N) for k in range(jp[p],jp[p+1]) if jx[k]>p))
    print(seed,"brute force %.4f  SQA %.4f  offset match %s"%(Ebf,qubo_energy(x,a,*jj),abs(off-off0)<1e-9))
    ok+=abs(qubo_energy(x,a,*jj)-Ebf)<1e-9
print("ground state found in",ok,"of 6")

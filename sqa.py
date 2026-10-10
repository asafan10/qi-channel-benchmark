"""
sqa.py -- sparse path-integral simulated quantum annealing (Martonak, Santoro, Tosatti 2002) on the QUBO's Ising form.

Ising form (Eq. 11 of the paper): E(s) = E0 + sum_p h_p s_p + sum_{p<q} Jc_pq s_p s_q,  Jc = J/4, h = a/2 + J.1/4.
Effective classical Hamiltonian on P Trotter slices at temperature T:
    H_eff = sum_k E(s^k) - Jperp(G) sum_k sum_p s_p^k s_p^{k+1},   Jperp = -(P T / 2) ln tanh(G / (P T)),
sampled with single-spin Metropolis at inverse temperature 1/(P T); G decreases linearly from G0 to G1.
The correct flip energy is used (fault 3 of the audited framework is not present):
    dE_cl = -2 s_p^k l_p^k,  l_p^k = h_p + sum_q Jc_pq s_q^k  (note: E(s) uses +h.s, so flipping s_p changes E by -2 s_p (h_p + sum Jc s));
    dE_q  = +2 Jperp s_p^k (s_p^{k-1} + s_p^{k+1}).
Returns the slice configuration with the lowest classical energy seen at the end of any sweep.
"""
import math
import numpy as np
from numba import njit


def ising_form(a, J):
    Jc = (J * 0.25).tocsr()
    h = a / 2 + np.asarray(J.sum(1)).ravel() / 4
    return h, Jc.indptr.astype(np.int64), Jc.indices.astype(np.int64), Jc.data.astype(np.float64)


@njit(cache=True)
def _energy(s, h, jp, jx, jd):
    e = 0.0
    for p in range(len(s)):
        e += h[p] * s[p]
        for k in range(jp[p], jp[p + 1]):
            q = jx[k]
            if q > p:
                e += jd[k] * s[p] * s[q]
    return e


@njit(cache=True)
def sqa(S0, h, jp, jx, jd, T, G0, G1, sweeps, seed):
    np.random.seed(seed)
    P, N = S0.shape; S = S0.copy()
    L = np.zeros((P, N))                      # local fields l_p^k
    for k in range(P):
        for p in range(N):
            v = h[p]
            for m in range(jp[p], jp[p + 1]):
                v += jd[m] * S[k, jx[m]]
            L[k, p] = v
    PT = P * T; beta = 1.0 / PT
    best = 1e300; bs = S[0].copy()
    for sw in range(sweeps):
        G = G0 + (G1 - G0) * sw / max(sweeps - 1, 1)
        Jp = -0.5 * PT * math.log(math.tanh(G / PT))
        for k in range(P):
            km = (k - 1) % P; kp = (k + 1) % P
            for _ in range(N):
                p = np.random.randint(N)
                sp = S[k, p]
                d = -2.0 * sp * L[k, p] + 2.0 * Jp * sp * (S[km, p] + S[kp, p])
                if d <= 0.0 or np.random.random() < math.exp(-beta * d):
                    S[k, p] = -sp
                    ds = -2.0 * sp
                    for m in range(jp[p], jp[p + 1]):
                        L[k, jx[m]] += jd[m] * ds
        for k in range(P):
            e = _energy(S[k], h, jp, jx, jd)
            if e < best:
                best = e; bs = S[k].copy()
    return bs, best


def sqa_solve(a, J, P, T, G0, G1, sweeps, rng, seed):
    h, jp, jx, jd = ising_form(a, J)
    S0 = rng.choice(np.array([-1.0, 1.0]), size=(P, len(a)))
    s, e = sqa(S0, h, jp, jx, jd, T, G0, G1, sweeps, seed)
    return ((s + 1) / 2).astype(np.int64), e

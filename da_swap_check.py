import numpy as np, math, itertools
rng = np.random.default_rng(1)
N = 12
A = rng.normal(0, 1, (N, N)); Q = (A + A.T) / 2                     # symmetric QUBO, E = x^T Q x
X = np.array(list(itertools.product([0, 1], repeat=N)))
E_all = np.einsum("ki,ij,kj->k", X, Q, X)
temps = np.logspace(math.log10(0.05), math.log10(3.0), 6)           # ascending, as in DA
exact = [float((E_all * np.exp(-(E_all - E_all.min()) / T)).sum() / np.exp(-(E_all - E_all.min()) / T).sum()) for T in temps]

def run(sign, sweeps=6000, seed=0):
    """sign=+1: criterion exactly as coded in DA (exp(dB*(E[r+1]-E[r]))); sign=-1: standard PT criterion."""
    rg = np.random.default_rng(seed); R = len(temps)
    x = rg.integers(0, 2, (R, N)); E = np.array([r @ Q @ r for r in x]); acc = np.zeros(R); cnt = 0; sw = 0; tot = 0
    for s in range(sweeps):
        for r in range(R):
            T = temps[r]
            for _ in range(N):
                i = rg.integers(0, N); xi = x[r, i]
                col = Q[:, i] + Q[i, :]; dE = (1 - 2 * xi) * (col @ x[r] - Q[i, i] * (2 * xi - 1))
                if dE < 0 or rg.random() < math.exp(max(-500.0, -dE / T)):
                    x[r, i] ^= 1; E[r] += dE
        if s % 3 == 0:
            for r in range(R - 1):
                dB = 1 / temps[r] - 1 / temps[r + 1]; dE = E[r + 1] - E[r]
                p = min(1.0, math.exp(max(-500, min(500.0, sign * dB * dE)))); tot += 1
                if rg.random() < p:
                    x[[r, r + 1]] = x[[r + 1, r]]; E[[r, r + 1]] = E[[r + 1, r]]; sw += 1
        if s > sweeps // 5: acc += E; cnt += 1
    return acc / cnt, sw / tot

print("T          :", np.round(temps, 3))
print("exact <E>  :", np.round(exact, 3))
for sign, name in ((+1, "as coded (exp(+dB*(E[r+1]-E[r])))"), (-1, "standard PT      ")):
    m, sr = run(sign); print(f"{name}: <E>", np.round(m, 3), f"swap acceptance {sr:.2f}  max|<E>-exact| = {np.abs(m-np.array(exact)).max():.3f}")

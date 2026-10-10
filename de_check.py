import numpy as np
rng = np.random.default_rng(0)
def true_dE(Q, x, i):
    y = x.copy(); y[i] ^= 1
    return float(y @ Q @ y - x @ Q @ x)
def code_dE(Q, x, i):                       # exactly as in DA._dE and inlined in PA.solve
    col = Q[:, i] + Q[i, :]
    return float((1 - 2 * x[i]) * (col @ x - Q[i, i] * x[i]))
for name in ("upper-triangular Q", "symmetric Q"):
    err0 = err1 = 0.0
    for _ in range(2000):
        N = 8; A = rng.normal(0, 1, (N, N)); Q = np.triu(A) if name.startswith("upper") else (A + A.T) / 2
        Q[np.diag_indices(N)] = rng.normal(-3, 1, N)               # large negative diagonal as in the one-hot penalty
        x = rng.integers(0, 2, N); i = rng.integers(0, N)
        e = abs(code_dE(Q, x, i) - true_dE(Q, x, i))
        if x[i] == 0: err0 = max(err0, e)
        else: err1 = max(err1, e)
    print(f"{name}: max |code dE - true dE|  for 0->1 flips: {err0:.3f}   for 1->0 flips: {err1:.3e}")

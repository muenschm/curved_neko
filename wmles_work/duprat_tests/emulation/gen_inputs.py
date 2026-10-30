"""Random wall-node inputs covering ZPG, FPG, APG and reversed-stress cases."""
import numpy as np
rng = np.random.default_rng(42)
n = 3000
nrm = rng.normal(size=(n, 3)); nrm /= np.linalg.norm(nrm, axis=1)[:, None]
nrm[:200] = [0, 1, 0]
vel = rng.normal(size=(n, 3)) * rng.choice([0.01, 0.1, 1.0], size=(n, 1))
vel[200:220] = nrm[200:220] * 0.5          # purely wall-normal
nu = 10**rng.uniform(-5, -3.5, n)
rho = rng.uniform(0.8, 1.5, n)
h = 10**rng.uniform(-2.3, -0.7, n)
dpds = rng.choice([-1, 1], n) * 10**rng.uniform(-4, 1, n)
dpds[220:300] = 0.0
grad = rng.normal(size=(n, 3)) * 10**rng.uniform(-3, 2, (n, 1))
tau_guess = rng.normal(size=(n, 3)) * 1e-3
with open('inputs.txt', 'w') as f:
    f.write(f'{n}\n')
    for i in range(n):
        vals = [*vel[i], *nrm[i], nu[i], rho[i], h[i], dpds[i], *grad[i], *tau_guess[i]]
        f.write(' '.join(f'{v:.17e}' for v in vals) + '\n')

import sys, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import warnings; warnings.simplefilter('ignore')
from duprat_reference import u_star, solve
# tau along u-hat as a function of u, hill-like numbers
nu = 1/10595
for h, dp in [(0.03, 0.05), (0.03, 0.5), (0.1, 0.2), (0.01, 1.0)]:
    up = (nu*dp)**(1/3)
    usep = up*u_star(h*up/nu, 0.0, 1, 1)
    print(f'h {h} dp/ds {dp}: up {up:.4f} u_sep {usep:.4f}')
    for u in np.array([0.02, 0.1, 0.3, 0.6, 0.9, 0.99, 1.01, 1.1, 1.5, 3.0])*usep:
        ut, al, st, _ = solve(u, h, nu, dp)
        print(f'   u {u:.4f} tau.uhat {st*ut**2:+.3e} alpha {al:.3f}')
    # count sign changes of the residual on each branch
    for st in (1.0, -1.0):
        for u in [0.3*usep, 1.5*usep]:
            uts = np.geomspace(1e-6, 1.0, 400)
            r = []
            for ut in uts:
                utp = np.sqrt(ut**2+up**2); al = ut**2/(ut**2+up**2)
                r.append(st*(utp*u_star(h*utp/nu, al, 1, st)-u))
            r = np.array(r)
            print(f'   branch {st:+.0f} u/u_sep {u/usep:.1f}: sign changes {np.sum(np.diff(np.sign(r))!=0)}')

#!/usr/bin/env python3
"""Reference values for the unit tests of the Duprat wall model.

U*(y*) is integrated with adaptive quadrature from Eq. (5) of Duprat et al.
(2011), including the sign of the wall shear stress,
  dU*/dy* = (s_p (1 - alpha)^(3/2) y* + s_tau alpha) / (1 + nu_t / nu),
with the eddy viscosity of Eq. (6). u_tau solves u = u_tau_p U*(h u_tau_p / nu)
with brentq. The branch s_tau = -1 is taken under an adverse gradient when
the separation profile (alpha = 0) already exceeds the sampled velocity.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

KAPPA, BETA, A = 0.41, 0.78, 17.0


def u_star(ys, alpha, sp, st, kappa=KAPPA, beta=BETA, a=A):
    apg = max(1.0 - alpha, 0.0)**1.5
    damp = 1.0 + a*alpha**3

    def f(eta):
        nut = kappa*eta*(alpha + eta*apg)**beta*(1 - np.exp(-eta/damp))**2
        return (sp*apg*eta + st*alpha)/(1 + nut)
    if ys <= 0:
        return 0.0
    # integrate in s = ln(1 + y*) like the kernel, but adaptively
    g = lambda s: f(np.expm1(s))*np.exp(s)
    v, _ = quad(g, 0.0, np.log1p(ys), epsabs=0, epsrel=2e-14, limit=400)
    return v


def solve(u, h, nu, dpds, rho=1.0, kappa=KAPPA, beta=BETA, a=A):
    up = (nu*abs(dpds)/rho)**(1/3)
    sp = 1.0 if dpds >= 0 else -1.0
    st = 1.0
    if sp > 0 and up > 0 and up*u_star(h*up/nu, 0.0, sp, st, kappa, beta, a) > u:
        st = -1.0

    def res(ut):
        utp = np.sqrt(ut**2 + up**2)
        al = ut**2/(ut**2 + up**2)
        return st*(utp*u_star(h*utp/nu, al, sp, st, kappa, beta, a) - u)
    lo, hi = 1e-12, 1.0
    while res(hi) <= 0:
        hi *= 2
    ut = brentq(res, lo, hi, xtol=1e-17, rtol=8.9e-16, maxiter=500)
    alpha = ut**2/(ut**2 + up**2) if up > 0 else 1.0
    return ut, alpha, st, up


if __name__ == '__main__':
    np.set_printoptions(precision=17)
    print('# checks against the existing references')
    for args, ref in [((1.0, 0.05, 1e-4, 0.0), 7.07651223726671053e-02),
                      ((1.0, 0.05, 1e-4, 1.0), 6.51766461571228473e-02),
                      ((1.0, 0.05, 1e-4, -1.0), 1.41161919435572131e-01),
                      ((0.5, 0.02, 1e-4, 0.2), 5.12131551453846173e-02)]:
        ut, al, st, up = solve(*args)
        print(args, f'{ut:.17e}', f'rel {abs(ut - ref)/ref:.1e}', st)
    print(f'U*(100, 0.5, +1, +1) = {u_star(100, 0.5, 1, 1):.17e}'
          ' (ref 1.61004369720788532e+01)')
    print('# new: U* on the reversed branch')
    for ys, al in [(100.0, 0.5), (50.0, 0.05), (5.0, 0.3)]:
        print(f'U*({ys}, {al}, +1, -1) = {u_star(ys, al, 1, -1):.17e}')
    print('# new: reversed wall stress')
    for args in [(0.1, 0.05, 1e-4, 10.0, 1.0), (0.1, 0.05, 1e-4, 12.0, 1.2),
                 (0.5, 0.05, 1e-4, 2.0, 1.0)]:
        ut, al, st, up = solve(*args)
        usep = up*u_star(args[1]*up/args[2], 0.0, 1.0, 1.0)
        print(args, f'utau {ut:.17e} alpha {al:.17e} st {st:+.0f} '
              f'up {up:.6f} u_sep {usep:.6f}')

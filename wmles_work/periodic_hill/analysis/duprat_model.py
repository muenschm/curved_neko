#!/usr/bin/env python3
"""The wall model of Duprat et al. (2011) in numpy, for post-processing.

Same equations as Neko's kernel (src/wall_models/bcknd/cpu/duprat_cpu.f90,
branch feature/duprat-local-dpds):

  u_p = (nu |dp/ds| / rho)^(1/3),  u_tau_p = sqrt(u_tau^2 + u_p^2),
  alpha = u_tau^2 / u_tau_p^2,
  dU*/dy* = (s_p (1 - alpha)^(3/2) y* + s_tau alpha) / (1 + nu_t / nu),
  nu_t / nu = kappa y* (alpha + y* (1 - alpha)^(3/2))^beta
              (1 - exp(-y* / (1 + A alpha^3)))^2,

with s_p the sign of dp/ds and s_tau the sign of the wall shear stress, both
relative to the flow at the sampling point. U* is integrated with the same
composite Gauss rule in ln(1 + y*) as the kernel. u_tau solves
u = u_tau_p U*(h u_tau_p / nu).

Branch: under an adverse gradient, the separation profile (alpha = 0) gives
u_sep = u_p U*(h u_p / nu). If the sampled u is smaller, the near-wall flow
is reversed, s_tau = -1.

The kernel uses a safeguarded Newton method; here plain bisection is used,
which is slower but simple. The root is unique on each branch for the
parameters of the hill (see duprat_tests/branch_scan.py).
"""
import numpy as np

KAPPA, BETA, A = 0.41, 0.78, 17.0
N_SUB = 12
_x, _w = np.polynomial.legendre.leggauss(5)
GL_X, GL_W = 0.5*(_x + 1.0), 0.5*_w          # 5-point Gauss on [0, 1]


def u_star(ys, alpha, sp, st, kappa=KAPPA, beta=BETA, a=A):
    """Dimensionless velocity U*(y*), Eq. (5) integrated from 0 to ys."""
    if ys <= 0.0:
        return 0.0
    apg = max(1.0 - alpha, 0.0)**1.5
    damp = 1.0 + a*alpha**3
    ds = np.log1p(ys)/N_SUB
    s = (np.arange(N_SUB)[:, None] + GL_X[None, :])*ds   # (N_SUB, 5)
    eta = np.expm1(s)
    nut = kappa*eta*(alpha + eta*apg)**beta*(1.0 - np.exp(-eta/damp))**2
    # dU*/dy* times the Jacobian d(eta)/ds = 1 + eta
    return float(np.sum(GL_W[None, :]*ds*(sp*apg*eta + st*alpha)
                        / (1.0 + nut)*(1.0 + eta)))


def friction_velocity(u, h, nu, dpds, rho=1.0, kappa=KAPPA, beta=BETA, a=A):
    """u_tau, alpha and s_tau for the wall-parallel speed u > 0 at distance h.

    dpds is the pressure gradient along the flow, positive when adverse.
    The wall shear stress along the flow direction is -s_tau rho u_tau^2
    (it opposes the near-wall flow, which is reversed for s_tau = -1)."""
    up = (nu*abs(dpds)/rho)**(1.0/3.0)
    sp = 1.0 if dpds >= 0.0 else -1.0
    st = 1.0
    if sp > 0 and up > 0 and up*u_star(h*up/nu, 0.0, sp, 1.0, kappa, beta,
                                       a) > u:
        st = -1.0

    def f(ut):
        # grows with ut on both branches
        utp = np.hypot(ut, up)
        al = ut**2/(ut**2 + up**2)
        return st*(utp*u_star(h*utp/nu, al, sp, st, kappa, beta, a) - u)

    lo = 1e-10*max(u, up)
    if f(lo) >= 0.0:                       # vanishing wall shear stress
        return lo, (lo**2/(lo**2 + up**2) if up > 0 else 1.0), st
    hi = max(u, up)
    while f(hi) <= 0.0:
        hi *= 2.0
    for _ in range(100):                   # bisection in log(u_tau)
        mid = np.sqrt(lo*hi)
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
        if hi/lo - 1.0 < 1e-12:
            break
    ut = np.sqrt(lo*hi)
    alpha = ut**2/(ut**2 + up**2) if up > 0 else 1.0
    return ut, alpha, st


if __name__ == '__main__':
    # Check against the SciPy references of the unit tests
    # (tests/unit/duprat/test_duprat.pf)
    for args, ref, st_ref in [((1.0, 0.05, 1e-4, 0.0), 7.07651223726671053e-02, 1),
                              ((1.0, 0.05, 1e-4, 1.0), 6.51766461571228473e-02, 1),
                              ((1.0, 0.05, 1e-4, -1.0), 1.41161919435572131e-01, 1),
                              ((0.1, 0.05, 1e-4, 10.0), 2.16429977179458027e-01, -1),
                              ((0.24, 0.03, 1e-4, 0.5), 5.86134660668170913e-03, -1),
                              ((0.245, 0.03, 1e-4, 0.5), 5.72722076270732253e-03, 1)]:
        ut, al, st = friction_velocity(*args)
        print(f'u, h, nu, dp/ds = {args}: u_tau = {ut:.10e} '
              f'(ref {ref:.10e}, rel. diff {abs(ut - ref)/ref:.1e}), '
              f's_tau = {st:+.0f} (ref {st_ref:+d})')

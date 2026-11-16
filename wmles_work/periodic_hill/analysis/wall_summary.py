#!/usr/bin/env python3
"""Separation, reattachment and skin friction on the hill wall.

Usage:
  python3 -I wall_summary.py RUNDIR START [--from-time T] [--min-window 40]
                             [--re 10595] [--sampling 3] [--name NAME]
                             [--out FILE.csv]

RUNDIR holds fluid_stats2d* and, for Duprat runs, wall_stats2d* (see
README.md). START as for ph_stats.py, e.g. 150 or 150,489.53 after a
restart. Windows are averaged as in ph_stats.py.

At every node of the hill wall the script takes the mean velocity at the
sampling point (GLL row --sampling of the wall element, as in the case file)
and its component u_t along the wall (positive in +x). It reports:

1. Zero crossings of <u_t>: separation (sep) and reattachment (reatt). This
   is the method of the report: the reconstructed C_f there has the sign of
   the mean flow at the sampling point.
2. For Duprat runs, C_f reconstructed with the Duprat model (duprat_model.py)
   from <u_t>, h and <dp/ds>: C_f = 2 sign(u_t) s_tau u_tau^2 (rho = 1,
   U_b = 1). Its sign also accounts for reversed near-wall flow under an
   adverse gradient (s_tau = -1), so its zero crossings can differ from 1.
3. C_f from the statistics of |tau_w|, signed with <u_t>: 2 sign(u_t)
   <|tau_w|>. Note <|tau_w|> is the mean magnitude, which is larger than the
   magnitude of the mean where the stress changes sign in time.
4. The mean <|tau_w|> over the hill wall (integral over x divided by 9),
   the alpha spike ratio at element boundaries (about 1 without spikes),
   and the mean C_f of 2. and 3. over x, to compare with Fx_v of
   force_torque.
5. The main bubble per window (windows of at least --min-window): mean and
   standard error of separation and reattachment.

--out writes x, y, <u_t>, h, <dp/ds>, <alpha>, C_f (model), C_f (|tau_w|)
along the wall.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from duprat_model import friction_velocity  # noqa: E402
from ph_io import LX, combine, hill_height, stats_files, wall_fields  # noqa


def wall_nodes(geo, js):
    """Wall nodes of the hill: element, GLL index, x, y, sampling point,
    unit tangent t (along +x), and the sampling distance h (projection of
    the offset to the sampling point on the wall normal, as in Neko)."""
    X, Y, lx = geo['x'], geo['y'], geo['lx']
    rows = []
    for e in range(X.shape[0]):
        yw = np.array([hill_height(v) for v in X[e, 0, :]])
        if np.all(np.abs(Y[e, 0, :] - yw) < 1e-5):
            for i in range(lx):
                rows.append((X[e, 0, i], Y[e, 0, i], e, i))
    rows.sort()
    # keep one of the two copies of a node shared by two elements
    keep, last = [], None
    for r in rows:
        if last is None or abs(r[0] - last) > 1e-7:
            keep.append(r)
            last = r[0]
    xw = np.array([r[0] for r in keep])
    yw = np.array([r[1] for r in keep])
    el = np.array([r[2] for r in keep])
    ix = np.array([r[3] for r in keep])
    tx, ty = np.gradient(xw), np.gradient(yw)
    tn = np.hypot(tx, ty)
    tx, ty = tx/tn, ty/tn
    xs, ys = X[el, js, ix], Y[el, js, ix]
    h = (xs - xw)*(-ty) + (ys - yw)*tx      # inward normal (-t_y, t_x)
    return dict(x=xw, y=yw, e=el, i=ix, tx=tx, ty=ty, h=h)


def u_tangential(nodes, f, js):
    e, i = nodes['e'], nodes['i']
    return f['u'][e, js, i]*nodes['tx'] + f['v'][e, js, i]*nodes['ty']


def crossings(x, q):
    """x of the sign changes of q, linear in between; 'sep' for + to -."""
    out = []
    for j in np.where(np.diff(np.sign(q)) != 0)[0]:
        xc = x[j] - q[j]*(x[j + 1] - x[j])/(q[j + 1] - q[j])
        out.append((xc, 'sep' if q[j] > 0 else 'reatt'))
    return out


def main_bubble(cr):
    """First separation in 0 <= x < 4 and the next reattachment."""
    sep = [c for c, k in cr if k == 'sep' and 0.0 <= c < 4.0]
    if not sep:
        return np.nan, np.nan
    reatt = [c for c, k in cr if k == 'reatt' and c > sep[0]]
    return sep[0], (reatt[0] if reatt else np.nan)


def cf_model(nodes, ut, wall, nu):
    """C_f reconstructed with the Duprat model from the mean flow."""
    cf = np.zeros(len(ut))
    for k in range(len(ut)):
        e, i = nodes['e'][k], nodes['i'][k]
        u = abs(ut[k])
        if u < 1e-12:
            continue
        utau, _, st = friction_velocity(u, nodes['h'][k], nu,
                                        wall['dpds'][e, 0, i])
        cf[k] = 2.0*np.sign(ut[k])*st*utau**2
    return cf


def fmt(cr):
    return ', '.join(f'{x:.2f}({k})' for x, k in cr) or '-'


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('rundir')
    ap.add_argument('start', help='start time(s), comma separated')
    ap.add_argument('--from-time', type=float, default=None)
    ap.add_argument('--min-window', type=float, default=40.0)
    ap.add_argument('--re', type=float, default=10595.0)
    ap.add_argument('--sampling', type=int, default=3)
    ap.add_argument('--name', default=None)
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    start = [float(v) for v in args.start.split(',')]
    nu, js = 1.0/args.re, args.sampling
    name = args.name or os.path.basename(os.path.normpath(args.rundir))

    geo, mean, used, _ = combine(stats_files(args.rundir), start,
                                 args.from_time)
    nodes = wall_nodes(geo, js)
    x = nodes['x']
    ut = u_tangential(nodes, mean, js)
    print(f'{name}: {len(used)} windows, {used[0][0]:g}-{used[-1][1]:g}; '
          f'{len(x)} wall nodes, h = {nodes["h"].min():.4f}'
          f'...{nodes["h"].max():.4f}')
    cr_u = crossings(x, ut)
    print(f'  1. zero crossings of <u_t>:        {fmt(cr_u)}')

    try:
        wseries = stats_files(args.rundir, 'wall_stats2d')
    except SystemExit:
        wseries = None
    cf_m = cf_t = None
    if wseries:
        _, wmean, wused, _ = combine(wseries, start, args.from_time,
                                     names=wall_fields)
        cf_m = cf_model(nodes, ut, wmean, nu)
        e, i = nodes['e'], nodes['i']
        cf_t = 2.0*np.sign(ut)*wmean['tau'][e, 0, i]
        bub = (x > 0.5) & (x < 4.0)
        ww = (x > 8.0) & (x < 9.0)
        print(f'  2. zero crossings of C_f (model):  {fmt(crossings(x, cf_m))}'
              f'; C_f min (bubble) {cf_m[bub].min():+.4f}, '
              f'max (windward) {cf_m[ww].max():.4f} at x = '
              f'{x[ww][np.argmax(cf_m[ww])]:.2f}')
        print(f'  3. C_f from <|tau_w|>: min (bubble) {cf_t[bub].min():+.4f}'
              f', max (windward) {cf_t[ww].max():.4f}')
        tauw = np.trapezoid(wmean['tau'][e, 0, i], x)/LX
        # alpha spikes: element-boundary nodes against interior nodes
        a = wmean['alpha'][e, 0, i]
        edge = (nodes['i'] == 0) | (nodes['i'] == geo['lx'] - 1)
        dev = np.abs(a[1:-1] - (a[:-2] + (a[2:] - a[:-2])*(x[1:-1] - x[:-2])
                                / (x[2:] - x[:-2])))
        ratio = np.median(dev[edge[1:-1]])/np.median(dev[~edge[1:-1]])
        print(f'  4. mean <|tau_w|> hill {tauw:.5f}, '
              f'alpha spike ratio {ratio:.1f}')
        # x component of the mean skin friction, comparable with Fx_v of
        # force_torque (tau_w t_x ds = tau_w dx)
        print(f'     mean C_f (x component): from <|tau_w|> '
              f'{np.trapezoid(cf_t, x)/LX:.5f}, model from mean flow '
              f'{np.trapezoid(cf_m, x)/LX:.5f}')

    # 5. per window
    seps, reatts = [], []
    for a_, b_, f in used:
        if b_ - a_ < args.min_window:
            continue
        s_, r_ = main_bubble(crossings(x, u_tangential(nodes, f, js)))
        seps.append(s_)
        reatts.append(r_)
    if len(seps) >= 2:
        s_, r_ = np.array(seps), np.array(reatts)
        se = lambda v: np.nanstd(v, ddof=1)/np.sqrt(np.sum(~np.isnan(v)))
        print(f'  5. per window ({len(s_)} of >= {args.min_window:g}): '
              f'sep {np.nanmean(s_):.2f} +- {se(s_):.2f}, '
              f'reatt {np.nanmean(r_):.2f} +- {se(r_):.2f} '
              f'(range {np.nanmin(r_):.2f}...{np.nanmax(r_):.2f})')

    if args.out:
        cols = [x, nodes['y'], ut, nodes['h']]
        hdr = 'x,y,u_t,h'
        if wseries:
            e, i = nodes['e'], nodes['i']
            cols += [wmean['dpds'][e, 0, i], wmean['alpha'][e, 0, i], cf_m,
                     cf_t]
            hdr += ',dpds,alpha,cf_model,cf_tau'
        np.savetxt(args.out, np.column_stack(cols), delimiter=',',
                   header=hdr, comments='')


if __name__ == '__main__':
    main()

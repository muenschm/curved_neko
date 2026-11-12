#!/usr/bin/env python3
"""Mean profiles of the periodic hill from Neko's z-averaged statistics.

Usage:
  python3 -I ph_stats.py DIR START [--from-time T] [--x 0.05,0.5,...]
                         [--out FILE.csv]

DIR holds the 2D outputs of fluid_stats ("avg_direction": "z",
"set_of_stats": "basic"), named <prefix>.fNNNNN. Neko writes the average
over the time since the previous output (or since the start of sampling)
and then resets the statistics, so every file is one averaging window.
START is the time sampling started; after a restart, list the restart time
too, e.g. 150,483.3943, since a restarted run starts a new window there.
The script averages the windows (weighted by their length), optionally only
those starting at or after --from-time, and writes the mean velocity and the
Reynolds stresses along vertical lines at the given x/h.

Field layout of the 2D files: U = (<u>, <v>), P = <p>, T = <pp>,
S = (<uu>, <vv>, <ww>, <uv>, <uw>, <vw>, <w>).

The module also provides the helpers used by gpu_cpu_windows.py.
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wall_dpds_check import hill_height, read_fld  # noqa: E402,F401

LX = 9.0
LY = 3.035
XS_DEFAULT = [0.05, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
S_NAMES = ['uu', 'vv', 'ww', 'uv', 'uw', 'vw', 'w']


def gll(n):
    """The n Gauss-Lobatto-Legendre points on [-1, 1]."""
    if n == 2:
        return np.array([-1.0, 1.0])
    p = np.polynomial.legendre.Legendre.basis(n - 1).deriv()
    return np.concatenate([[-1.0], np.sort(p.roots().real), [1.0]])


class Lagrange:
    """Lagrange interpolation on given nodes (barycentric form)."""

    def __init__(self, nodes):
        self.x = np.asarray(nodes, float)
        d = self.x[:, None] - self.x[None, :]
        np.fill_diagonal(d, 1.0)
        self.w = 1.0/np.prod(d, axis=1)

    def basis(self, r):
        """Basis values, shape (len(r), n)."""
        r = np.atleast_1d(np.asarray(r, float))
        d = r[:, None] - self.x[None, :]
        hit = np.isclose(d, 0.0, atol=1e-14)
        d[hit] = 1.0
        t = self.w/d
        phi = t/np.sum(t, axis=1, keepdims=True)
        rows = np.any(hit, axis=1)
        phi[rows] = hit[rows].astype(float)
        return phi

    def dbasis(self, r, h=1e-7):
        """Derivative of the basis, by central differences."""
        return (self.basis(np.asarray(r) + h)
                - self.basis(np.asarray(r) - h))/(2*h)


def stats_files(d, prefix='fluid_stats2d'):
    """The statistics files in d, as a list of series, each sorted by
    number. A series is a set of files with the same name before .f."""
    fs = glob.glob(os.path.join(d, prefix + '*.f[0-9]*'))
    if not fs:
        raise SystemExit(f'no {prefix}*.f* in {d}')
    series = {}
    for f in fs:
        m = re.match(r'(.*)\.f(\d+)$', os.path.basename(f))
        series.setdefault(m.group(1), []).append((int(m.group(2)), f))
    return [[f for _, f in sorted(v)] for _, v in sorted(series.items())]


def fields(r):
    """The statistics of one file as a dict of (nel, ly, lx) arrays."""
    f = {'u': r['U'][:, 0], 'v': r['U'][:, 1], 'p': r['P'], 'pp': r['T']}
    for k, name in enumerate(S_NAMES):
        f[name] = r['S'][:, k]
    return f


def combine(series, start, from_time=None):
    """Read all windows and average them.

    Returns (geo, mean, windows, info): geo has the coordinates 'x', 'y'
    (nel, ly, lx) and 'lx'; mean is the length-weighted average of the
    windows used; windows is a list of (t_begin, t_end, fields); info lists
    the windows skipped and why."""
    start = sorted(start)
    geo, windows, info = None, [], []
    for files in series:
        prev = None
        for f in files:
            r = read_fld(f)
            if 'X' in r and geo is None:
                geo = {'x': r['X'][:, 0], 'y': r['X'][:, 1], 'lx': r['lx'],
                       'ids': r['ids']}
            b = r['time']
            begins = [s for s in start if s < b]
            if not begins:
                info.append(f'{os.path.basename(f)}: t = {b:g} before start')
                prev = b
                continue
            a = max(begins[-1], prev if prev is not None else -np.inf)
            prev = b
            if not any(np.any(r[k] != 0) for k in ('U', 'S') if k in r):
                info.append(f'{os.path.basename(f)}: empty')
                continue
            if b - a <= 1e-6*max(1.0, abs(b)):
                info.append(f'{os.path.basename(f)}: zero length')
                continue
            windows.append((a, b, fields(r)))
    if geo is None:
        raise SystemExit('no coordinates found in the statistics files')
    windows.sort(key=lambda w: w[0])
    used = [w for w in windows if from_time is None or w[0] >= from_time]
    if not used:
        raise SystemExit('no averaging window left')
    tot = sum(b - a for a, b, _ in used)
    mean = {k: sum((b - a)*f[k] for a, b, f in used)/tot
            for k in used[0][2]}
    return geo, mean, used, info


def reynolds(f):
    """Mean velocity, pressure and Reynolds stresses from the raw moments."""
    return {'u': f['u'], 'v': f['v'], 'w': f['w'], 'p': f['p'],
            'upup': f['uu'] - f['u']**2, 'vpvp': f['vv'] - f['v']**2,
            'wpwp': f['ww'] - f['w']**2, 'upvp': f['uv'] - f['u']*f['v'],
            'k': 0.5*(f['uu'] - f['u']**2 + f['vv'] - f['v']**2
                      + f['ww'] - f['w']**2)}


def locate(geo, lag, x, y, tol=1e-10, slack=1e-3):
    """Element and reference coordinates (e, r, s) of the points (x, y).

    For every element whose bounding box contains the point, (r, s) is found
    with Newton's method on x(r, s) = sum_ij l_i(r) l_j(s) x_ij. The element
    with the smallest overshoot max(|r|, |s|) - 1 is taken. A small overshoot
    (up to `slack`) is accepted and clamped: the coordinates in the files are
    single precision, and the analytic hill surface differs slightly from the
    polynomial one, so points on element edges or on the wall can land just
    outside every element. Points further outside get e = -1."""
    X, Y = geo['x'], geo['y']
    xmin, xmax = X.min(axis=(1, 2)), X.max(axis=(1, 2))
    ymin, ymax = Y.min(axis=(1, 2)), Y.max(axis=(1, 2))
    n = len(x)
    e_out, r_out, s_out = -np.ones(n, int), np.zeros(n), np.zeros(n)
    for k in range(n):
        pad = 1e-3
        cand = np.where((x[k] >= xmin - pad) & (x[k] <= xmax + pad)
                        & (y[k] >= ymin - pad) & (y[k] <= ymax + pad))[0]
        best, best_rs, best_over = -1, None, np.inf
        for e in cand:
            rs = np.zeros(2)
            for _ in range(30):
                pr, ps = lag.basis(rs[0])[0], lag.basis(rs[1])[0]
                dr, ds = lag.dbasis(rs[0])[0], lag.dbasis(rs[1])[0]
                fx = ps @ X[e] @ pr - x[k]
                fy = ps @ Y[e] @ pr - y[k]
                J = np.array([[ps @ X[e] @ dr, ds @ X[e] @ pr],
                              [ps @ Y[e] @ dr, ds @ Y[e] @ pr]])
                step = np.linalg.solve(J, [fx, fy])
                rs = np.clip(rs - step, -1.5, 1.5)
                if np.max(np.abs(step)) < tol:
                    break
            over = np.max(np.abs(rs)) - 1.0
            if over < best_over:
                best, best_rs, best_over = e, rs, over
            if over <= 0.0:
                break
        if best >= 0 and best_over <= slack:
            rs = np.clip(best_rs, -1.0, 1.0)
            e_out[k], r_out[k], s_out[k] = best, rs[0], rs[1]
    return e_out, r_out, s_out


def interp(field, lag, e, r, s):
    """Values of a (nel, ly, lx) field at located points; NaN outside."""
    out = np.full(len(e), np.nan)
    for k in np.where(e >= 0)[0]:
        out[k] = lag.basis(s[k])[0] @ field[e[k]] @ lag.basis(r[k])[0]
    return out


def profiles(geo, mean, xs, n=201):
    """Profiles of the mean quantities along x = const, from the hill
    surface to the top wall."""
    lag = Lagrange(gll(geo['lx']))
    m = reynolds(mean)
    rows = []
    for x in xs:
        yy = np.linspace(hill_height(x), LY, n)
        e, r, s = locate(geo, lag, np.full_like(yy, x), yy)
        vals = {q: interp(m[q], lag, e, r, s) for q in
                ('u', 'v', 'upup', 'vpvp', 'wpwp', 'upvp', 'k')}
        for j in range(n):
            rows.append([x, yy[j]] + [vals[q][j] for q in
                                      ('u', 'v', 'upup', 'vpvp', 'wpwp',
                                       'upvp', 'k')])
    return np.array(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('dir')
    ap.add_argument('start', help='start time(s), comma separated')
    ap.add_argument('--from-time', type=float, default=None)
    ap.add_argument('--x', default=','.join(str(v) for v in XS_DEFAULT))
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    start = [float(v) for v in args.start.split(',')]
    geo, mean, used, info = combine(stats_files(args.dir), start,
                                    args.from_time)
    print(f'{args.dir}: {len(used)} windows, '
          + ', '.join(f'{a:g}-{b:g}' for a, b, _ in used))
    for line in info:
        print('  skipped', line)
    xs = [float(v) for v in args.x.split(',')]
    prof = profiles(geo, mean, xs)
    hdr = 'x,y,u,v,upup,vpvp,wpwp,upvp,k'
    for x in xs:
        p = prof[prof[:, 0] == x]
        print(f'  x/h = {x:4.2f}: min <u> {np.nanmin(p[:, 2]):+.4f}, '
              f'max u\'u\' {np.nanmax(p[:, 4]):.4f}, '
              f'max -u\'v\' {np.nanmax(-p[:, 7]):.4f}, '
              f'max v\'v\' {np.nanmax(p[:, 5]):.4f}, '
              f'points outside {int(np.sum(np.isnan(p[:, 2])))}')
    if args.out:
        np.savetxt(args.out, prof, delimiter=',', header=hdr, comments='')


if __name__ == '__main__':
    main()

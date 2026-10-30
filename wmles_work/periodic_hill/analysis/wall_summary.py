#!/usr/bin/env python3
"""Separation, wall shear stress and alpha spikes on the hill wall.

Usage:
  python3 -I wall_summary.py RUNDIR [RUNDIR ...]

RUNDIR holds fluid_stats2d* and wall_stats2d* as described in
wall_dpds_check.py. For each run the script prints:

- the x/h where the mean wall-parallel velocity at the sampling row (GLL
  index 3) changes sign, marked (sep) or (reatt);
- the mean <|tau_w|> over the hill wall, integrated in x and divided by the
  domain length 9;
- the alpha spike ratio: the median deviation of <alpha> from linear
  interpolation of its neighbours, at element-boundary nodes divided by the
  same at interior nodes (about 1 without spikes).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wall_dpds_check import hill_height, last_file, read_fld  # noqa: E402

JS = 3          # sampling row


def summary(d):
    fp, f0 = last_file(d, 'fluid_stats2d')
    wp, _ = last_file(d, 'wall_stats2d')
    fl, wl, X = read_fld(fp), read_fld(wp), read_fld(f0)['X']
    xs, ys = X[:, 0], X[:, 1]
    lx = fl['lx']
    bottom = [e for e in range(len(xs))
              if np.all(np.abs(ys[e, 0, :] - np.array(
                  [hill_height(v) for v in xs[e, 0, :]])) < 1e-6)]
    R = sorted((xs[e, 0, i], xs[e, JS, i], ys[e, JS, i], fl['U'][e, 0, JS, i],
                fl['U'][e, 1, JS, i], wl['S'][e, 0, 0, i], wl['U'][e, 0, 0, i],
                i in (0, lx - 1))
               for e in bottom for i in range(lx))
    R = np.array(R)
    _, k = np.unique(np.round(R[:, 0], 7), return_index=True)
    R = R[k]
    tx, ty = np.gradient(R[:, 1]), np.gradient(R[:, 2])
    tn = np.hypot(tx, ty)
    ut = (R[:, 3]*tx + R[:, 4]*ty)/tn
    sc = np.where(np.diff(np.sign(ut)) != 0)[0]
    cr = [f'{R[j, 0] - ut[j]*(R[j + 1, 0] - R[j, 0])/(ut[j + 1] - ut[j]):.2f}'
          f'{"(sep)" if ut[j] > 0 else "(reatt)"}' for j in sc]
    x, edge, a = R[:, 0], R[:, 7] > 0, R[:, 5]
    dev = np.abs(a[1:-1] - (a[:-2] + (a[2:] - a[:-2])*(x[1:-1] - x[:-2])
                            / (x[2:] - x[:-2])))
    e = edge[1:-1]
    tauw = np.trapezoid(R[:, 6], x)/9.0
    return (f'{os.path.basename(os.path.normpath(d)):8s} sign changes of '
            f'<u_t>: {", ".join(cr)} | mean <|tau_w|> hill {tauw:.5f} | '
            f'alpha spike ratio {np.median(dev[e])/np.median(dev[~e]):.1f}')


if __name__ == '__main__':
    for d in sys.argv[1:]:
        print(summary(d))

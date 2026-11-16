#!/usr/bin/env python3
"""Check the pressure gradient seen by the Duprat wall model on the hill.

Usage:
  python3 -I wall_dpds_check.py RUNDIR START [--from-time T] [--sampling 3]

RUNDIR holds fluid_stats2d* (fluid_stats, "avg_direction": "z") and
wall_stats2d* (user_stats of tau, duprat_dpds and duprat_alpha). START as for
ph_stats.py; the windows are averaged as there. The script reports, on the
hill wall:

- element boundaries: the largest difference between the two copies of a
  wall node shared by two neighbouring elements, for <dp/ds> and <alpha>;
- per region: the mean |<dp/ds>| used by the model, the mean |d<p>/ds| of
  the mean pressure along the row of sampling points (signed with the mean
  wall-parallel velocity there, as in the model), and their ratio.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ph_io import combine, hill_height, stats_files, wall_fields  # noqa: E402

REGIONS = [('crest', 0.0, 0.4), ('bubble', 0.5, 3.0),
           ('recovery', 4.0, 7.0), ('windward', 8.3, 8.8)]


def check(d, start, from_time, sampling):
    geo, fl, used, _ = combine(stats_files(d), start, from_time)
    _, wl, _, _ = combine(stats_files(d, 'wall_stats2d'), start, from_time,
                          names=wall_fields)
    ny = geo['x'].shape[1]
    xs, ys = geo['x'], geo['y']
    # bottom-wall elements: row j = 0 lies on the hill surface
    bottom = [e for e in range(len(xs))
              if np.all(np.abs(ys[e, 0, :] - np.array(
                  [hill_height(v) for v in xs[e, 0, :]])) < 1e-5)]
    js = sampling                        # sampling row, as in the model
    if js >= ny:
        raise SystemExit('sampling index outside the element')
    rows = []
    for e in bottom:
        for i in range(geo['lx']):
            xw, yw = xs[e, 0, i], ys[e, 0, i]
            xp, yp = xs[e, js, i], ys[e, js, i]
            up, vp = fl['u'][e, js, i], fl['v'][e, js, i]
            rows.append((xw, yw, xp, yp, up, vp, fl['p'][e, js, i],
                         wl['dpds'][e, 0, i], wl['alpha'][e, 0, i], e, i))
    r = np.array(rows)
    order = np.lexsort((r[:, 10], r[:, 0]))
    r = r[order]
    # element boundaries: equal wall positions from two elements
    same = np.where(np.abs(np.diff(r[:, 0])) < 1e-9)[0]
    jump_dpds = np.abs(r[same, 7] - r[same + 1, 7])
    jump_alpha = np.abs(r[same, 8] - r[same + 1, 8])
    scale = np.mean(np.abs(r[:, 7]))
    # unique points along the wall for the gradient of the mean pressure
    keep = np.ones(len(r), bool)
    keep[same + 1] = False
    u = r[keep]
    t = np.gradient(u[:, 0]), np.gradient(u[:, 1])
    tn = np.hypot(*t)
    tx, ty = t[0]/tn, t[1]/tn
    ut = u[:, 4]*tx + u[:, 5]*ty
    s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(u[:, 2]),
                                                 np.diff(u[:, 3])))])
    dpds_mean = np.gradient(u[:, 6], s)*np.sign(ut)
    lines = [f'{d}: {len(used)} windows, {used[0][0]:g}-{used[-1][1]:g}, '
             f'{len(bottom)} hill-wall elements, '
             f'{len(same)} element-boundary nodes',
             f'  element boundaries: max |jump <dp/ds>| = '
             f'{jump_dpds.max():.3e} (mean |<dp/ds>| = {scale:.3e}), '
             f'max |jump <alpha>| = {jump_alpha.max():.3e}']
    for name, a, b_ in REGIONS:
        m = (u[:, 0] >= a) & (u[:, 0] < b_)
        mod = np.mean(np.abs(u[m, 7]))
        mp = np.mean(np.abs(dpds_mean[m]))
        lines.append(f'  {name:9s} x/h {a:3.1f}-{b_:3.1f}: model '
                     f'<|<dp/ds>|> = {mod:.4f}, mean pressure <|d<p>/ds|> '
                     f'= {mp:.4f}, ratio model/mean = {mod/mp:.2f}, '
                     f'<alpha> = {np.mean(u[m, 8]):.3f}')
    return '\n'.join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('rundir')
    ap.add_argument('start', help='start time(s), comma separated')
    ap.add_argument('--from-time', type=float, default=None)
    ap.add_argument('--sampling', type=int, default=3)
    args = ap.parse_args()
    start = [float(v) for v in args.start.split(',')]
    print(check(args.rundir, start, args.from_time, args.sampling))


if __name__ == '__main__':
    main()

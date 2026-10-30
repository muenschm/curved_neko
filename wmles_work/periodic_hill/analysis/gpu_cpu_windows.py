#!/usr/bin/env python3
"""Compare two periodic-hill runs window by window, e.g. GPU against CPU.

Usage:
  python3 gpu_cpu_windows.py --run NAME DIR START --run NAME DIR START
                             [--min-window 40] [--from-time T] [--out FILE]

DIR and START are as for ph_stats.py (START may list one start time per
output series, e.g. 150,483.3943 after a restart). For every output window
of at least --min-window time units, the script takes the peaks of u'u',
-u'v' and v'v' and the minimum of <u> along the profiles at several x/h.
It reports the mean over the windows, its standard error (the scatter of
the window values divided by the square root of their number), and the
difference between the runs in percent and in units of the combined
standard error. The windows of one run are treated as independent samples.
"""
import argparse
import os

import numpy as np

from ph_stats import (LY, Lagrange, combine, gll, hill_height, interp,
                      locate, reynolds, stats_files)

XS = [0.5, 1.0, 2.0, 3.0, 4.0, 6.0]
QUANTITIES = [('max u\'u\'', 'upup', np.max, 1.0),
              ('max -u\'v\'', 'upvp', np.max, -1.0),
              ('max v\'v\'', 'vpvp', np.max, 1.0),
              ('min <u>', 'u', np.min, 1.0)]


def window_values(d, start, from_time, min_window):
    geo, _, windows, _ = combine(stats_files(d), start, from_time)
    lag = Lagrange(gll(geo['lx']))
    pts = {}
    for x in XS:
        yy = np.linspace(hill_height(x), LY, 201)
        pts[x] = locate(geo, lag, np.full_like(yy, x), yy)
    rows, used = [], []
    for a, b, f in windows:
        if b - a < min_window:
            continue
        m = reynolds(f)
        row = []
        for x in XS:
            e, r, s = pts[x]
            for _, q, op, sign in QUANTITIES:
                row.append(op(sign*interp(m[q], lag, e, r, s)))
        rows.append(row)
        used.append((a, b))
    return np.array(rows), used


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--run', nargs=3, action='append', required=True,
                    metavar=('NAME', 'DIR', 'START'))
    ap.add_argument('--min-window', type=float, default=40.0)
    ap.add_argument('--from-time', type=float, default=None)
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    if len(args.run) != 2:
        raise SystemExit('give exactly two runs')
    res = []
    for name, d, start in args.run:
        st = [float(v) for v in start.split(',')]
        vals, used = window_values(d, st, args.from_time, args.min_window)
        res.append((name, vals, used))
    (na, va, ua), (nb, vb, ub) = res
    lines = []
    for name, vals, used in res:
        lines.append(f'{name}: {len(used)} windows, ' + ', '.join(
            f'{a:g}-{b:g}' for a, b in used))
    lines.append(f'{"quantity":<19} {na:<20} {nb:<20} {"(A-B)/B":>8} '
                 f'{"(A-B)/sigma":>11}')
    k, higher, nturb = 0, 0, 0
    for x in XS:
        for label, q, _, _ in QUANTITIES:
            a, b = va[:, k], vb[:, k]
            ma, mb = a.mean(), b.mean()
            sa = a.std(ddof=1)/np.sqrt(len(a))
            sb = b.std(ddof=1)/np.sqrt(len(b))
            z = (ma - mb)/np.hypot(sa, sb)
            lines.append(f'x/h={x:<4} {label:<9} {ma:8.5f}+-{sa:.5f}  '
                         f'{mb:8.5f}+-{sb:.5f}  {100*(ma - mb)/abs(mb):+7.1f}%'
                         f'  {z:+6.1f}')
            if q != 'u':
                nturb += 1
                higher += ma > mb
            k += 1
    lines.append(f'turbulence peaks higher in {na}: {higher} of {nturb}')
    text = '\n'.join(lines)
    print(text)
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, 'w') as f:
            f.write(text + '\n')


if __name__ == '__main__':
    main()

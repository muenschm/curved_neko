#!/usr/bin/env python3
"""Time-averaged forces and forcing of a periodic-hill run.

Usage:
  python3 -I ph_forces.py DIR [--from-time 150] [--window 50]

DIR holds forces_zone3.csv (hill wall), forces_zone4.csv (top wall) and
forcing.csv, written by ../forces_from_log.py from the Neko logs. The
force_torque "scale" in the case file is 2 / (rho U_b^2 L_x L_z), so the
forces are coefficients: Fx_v is the mean skin-friction coefficient of a
wall and Fx_p its pressure-drag coefficient.

Averages are time-weighted (trapezoidal rule, the time step varies) over
t >= --from-time. The uncertainty is the standard error of the means over
windows of --window time units; these windows are not fully independent, so
it is a lower bound.

Balance: with the flow rate held constant, the body force f (per unit
volume) balances the x force on both walls on average,
<F_x,hill + F_x,top> = 2 H_mean <f>, with the mean channel height H_mean =
L_y - (1/L_x) int hill dx.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ph_io import LX, LY, hill_height  # noqa: E402


def load(path):
    d = np.genfromtxt(path, delimiter=',', names=True)
    _, k = np.unique(d['time'], return_index=True)      # drop duplicates
    return d[k]


def tavg(t, q):
    """Time-weighted average over the samples."""
    return np.trapezoid(q, t)/(t[-1] - t[0])


def mean_and_error(t, q, t0, window):
    m = t >= t0
    t, q = t[m], q[m]
    total = tavg(t, q)
    edges = list(np.arange(t0, t[-1] + 1e-9, window))
    if t[-1] - edges[-1] > 0.5*window:     # last window, nearly complete
        edges.append(t[-1])
    wm = []
    for a, b in zip(edges[:-1], edges[1:]):
        k = (t >= a) & (t <= b)
        if np.sum(k) > 2:
            wm.append(tavg(t[k], q[k]))
    err = np.std(wm, ddof=1)/np.sqrt(len(wm)) if len(wm) > 1 else np.nan
    return total, err, len(wm)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('dir')
    ap.add_argument('--from-time', type=float, default=150.0)
    ap.add_argument('--window', type=float, default=50.0)
    args = ap.parse_args()
    z3 = load(os.path.join(args.dir, 'forces_zone3.csv'))
    z4 = load(os.path.join(args.dir, 'forces_zone4.csv'))
    fo = load(os.path.join(args.dir, 'forcing.csv'))
    xs = np.linspace(0.0, LX, 9001)
    h_mean = LY - np.trapezoid([hill_height(x) for x in xs], xs)/LX

    t0, w = args.from_time, args.window
    rows = [('hill wall, pressure drag Fx_p', z3, 'Fx_p'),
            ('hill wall, skin friction Fx_v', z3, 'Fx_v'),
            ('top wall, skin friction Fx_v', z4, 'Fx_v')]
    print(f'{args.dir}: t >= {t0:g}, windows of {w:g}, '
          f'H_mean = {h_mean:.4f}')
    res = {}
    for label, d, col in rows:
        m, e, n = mean_and_error(d['time'], d[col], t0, w)
        res[label] = m
        print(f'  {label:32s} {m:.5f} +- {e:.5f}  ({n} windows)')
    mf, ef, n = mean_and_error(fo['time'], fo['f'], t0, w)
    print(f'  {"forcing <f>":32s} {mf:.5f} +- {ef:.5f}  ({n} windows)')
    fx3, _, _ = mean_and_error(z3['time'], z3['Fx'], t0, w)
    fx4, _, _ = mean_and_error(z4['time'], z4['Fx'], t0, w)
    print(f'  {"balance <Fx3 + Fx4>/(2 H f)":32s} '
          f'{(fx3 + fx4)/(2*h_mean*mf):.3f}')


if __name__ == '__main__':
    main()

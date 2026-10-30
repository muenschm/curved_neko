#!/usr/bin/env python3
"""Check the pressure gradient seen by the Duprat wall model on the hill.

Usage:
  python3 wall_dpds_check.py RUNDIR [RUNDIR ...] [--sampling 3]

RUNDIR holds the z-averaged statistics of one averaging window:
fluid_stats2d0.f00000 (fluid_stats, "avg_direction": "z") and
wall_stats2d0.f00000 (user_stats of tau, duprat_dpds and duprat_alpha). The
last file of each series that holds data is used. For each run the script
reports, on the hill wall:

- element boundaries: the largest difference between the two copies of a
  wall node shared by two neighbouring elements, for <dp/ds> and <alpha>;
- per region: the mean |<dp/ds>| used by the model, the mean |d<p>/ds| of
  the mean pressure along the row of sampling points (signed with the mean
  wall-parallel velocity there, as in the model), and their ratio.
"""
import argparse
import glob
import os
import struct

import numpy as np

LX_DOMAIN = 9.0
REGIONS = [('crest', 0.0, 0.4), ('bubble', 0.5, 3.0),
           ('recovery', 4.0, 7.0), ('windward', 8.3, 8.8)]


def hill_step(x, w=1.929, h=1.0):
    """Height of one half of the hill (Breuer et al.), as in the user file."""
    xs = x/w
    if xs <= 0:
        return h
    if xs <= 9/54:
        return h*min(1.0, 1 + 7.05575248e-1*xs**2 - 1.1947737203e1*xs**3)
    if xs <= 14/54:
        return h*(0.895484248 + 1.881283544*xs - 10.582126017*xs**2
                  + 10.627665327*xs**3)
    if xs <= 20/54:
        return h*(0.92128609 + 1.582719366*xs - 9.430521329*xs**2
                  + 9.147030728*xs**3)
    if xs <= 30/54:
        return h*(1.445155365 - 2.660621763*xs + 2.026499719*xs**2
                  - 1.164288215*xs**3)
    if xs <= 40/54:
        return h*(0.640164762 + 1.6863274926*xs - 5.798008941*xs**2
                  + 3.530416981*xs**3)
    if xs <= 1:
        return h*(2.013932568 - 3.877432121*xs + 1.713066537*xs**2
                  + 0.150433015*xs**3)
    return 0.0


def hill_height(x):
    x = x % LX_DOMAIN
    return hill_step(x) + hill_step(LX_DOMAIN - x)


def read_fld(path):
    """Read a Nek5000-format field file (2D or 3D, single or double)."""
    b = open(path, 'rb').read()
    hdr = b[:132].decode('ascii', 'replace').split()
    wdsz, lx, ly, lz = (int(v) for v in hdr[1:5])
    nel, time = int(hdr[5]), float(hdr[7])
    rdcode = hdr[-1]
    off = 132
    tag = struct.unpack('<f', b[off:off + 4])[0]
    endian = '<' if abs(tag - 6.54321) < 1e-5 else '>'
    off += 4
    ids = np.frombuffer(b, dtype=endian + 'i4', count=nel, offset=off)
    off += 4*nel
    dt = np.dtype(endian + ('f4' if wdsz == 4 else 'f8'))
    npt = lx*ly*lz
    ndim = 2 if lz == 1 else 3
    out = {'time': time, 'lx': lx, 'ly': ly, 'ids': ids.copy()}
    i = 0
    while i < len(rdcode):
        c = rdcode[i]
        if c in 'XU':
            a = np.frombuffer(b, dtype=dt, count=nel*ndim*npt, offset=off)
            out[c] = a.reshape(nel, ndim, ly, lx).astype(float)
            off += a.nbytes
            i += 1
        elif c in 'PT':
            a = np.frombuffer(b, dtype=dt, count=nel*npt, offset=off)
            out[c] = a.reshape(nel, ly, lx).astype(float)
            off += a.nbytes
            i += 1
        elif c == 'S':
            ns = int(rdcode[i + 1:i + 3])
            a = np.frombuffer(b, dtype=dt, count=nel*ns*npt, offset=off)
            out['S'] = a.reshape(nel, ns, ly, lx).astype(float)
            off += a.nbytes
            i += 3
        else:
            raise ValueError(f'{path}: unknown field code {c} in {rdcode}')
    return out


def last_file(d, prefix):
    """The last file of a series that holds data, and the first file, which
    holds the coordinates. A run may end with an empty output."""
    fs = sorted(glob.glob(os.path.join(d, prefix + '*.f[0-9]*')))
    if not fs:
        raise SystemExit(f'no {prefix}*.f* in {d}')
    for f in reversed(fs):
        r = read_fld(f)
        if any(np.any(v != 0) for k, v in r.items()
               if k in ('U', 'P', 'S')):
            return f, fs[0]
    raise SystemExit(f'all {prefix}*.f* in {d} are empty')


def check(d, sampling):
    fpath, ffirst = last_file(d, 'fluid_stats2d')
    wpath, wfirst = last_file(d, 'wall_stats2d')
    fl, wl = read_fld(fpath), read_fld(wpath)
    X = read_fld(ffirst)['X']            # coordinates are in the first file
    ny = fl['ly']
    xs, ys = X[:, 0], X[:, 1]
    # bottom-wall elements: row j = 0 lies on the hill surface
    bottom = [e for e in range(len(xs))
              if np.all(np.abs(ys[e, 0, :] - np.array(
                  [hill_height(v) for v in xs[e, 0, :]])) < 1e-6)]
    js = sampling                        # sampling row, as in the model
    if js >= ny:
        raise SystemExit('sampling index outside the element')
    rows = []
    for e in bottom:
        for i in range(fl['lx']):
            xw, yw = xs[e, 0, i], ys[e, 0, i]
            xp, yp = xs[e, js, i], ys[e, js, i]
            up, vp = fl['U'][e, 0, js, i], fl['U'][e, 1, js, i]
            rows.append((xw, yw, xp, yp, up, vp, fl['P'][e, js, i],
                         wl['U'][e, 1, 0, i], wl['S'][e, 0, 0, i], e, i))
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
    lines = [f'{d}: t = {fl["time"]:g}, {len(bottom)} hill-wall elements, '
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
    ap.add_argument('rundirs', nargs='+')
    ap.add_argument('--sampling', type=int, default=3)
    args = ap.parse_args()
    for d in args.rundirs:
        print(check(d, args.sampling))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Reading the periodic-hill statistics: geometry, file reader, windows.

Shared by the other scripts in this directory:

- `hill_height(x)`: the hill of Breuer et al., as in periodic_hill.f90;
- `read_fld(path)`: reads a Nek5000-format field file (2D or 3D, single or
  double precision);
- `stats_files(d, prefix)`: the statistics files in d, grouped into series
  (a restarted run starts a new series, e.g. fluid_stats2d00 and
  fluid_stats2d60);
- `fields(r)`, `wall_fields(r)`: the named fields of a fluid_stats or a
  wall_stats file;
- `combine(series, start, from_time)`: the averaging windows and their
  length-weighted mean.
"""
import glob
import os
import re
import struct

import numpy as np

LX = 9.0          # domain length
LY = 3.035        # domain height (top wall)
S_NAMES = ['uu', 'vv', 'ww', 'uv', 'uw', 'vw', 'w']


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
    x = x % LX
    return hill_step(x) + hill_step(LX - x)


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
            # each scalar is a separate block over all elements; out['S'] is
            # indexed [element, scalar, j, i]
            ns = int(rdcode[i + 1:i + 3])
            a = np.frombuffer(b, dtype=dt, count=nel*ns*npt, offset=off)
            out['S'] = a.reshape(ns, nel, ly, lx).transpose(1, 0, 2, 3) \
                .astype(float)
            off += a.nbytes
            i += 3
        else:
            raise ValueError(f'{path}: unknown field code {c} in {rdcode}')
    return out


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
    """The statistics of one fluid_stats file as a dict of (nel, ly, lx)
    arrays."""
    f = {'u': r['U'][:, 0], 'v': r['U'][:, 1], 'p': r['P'], 'pp': r['T']}
    for k, name in enumerate(S_NAMES):
        f[name] = r['S'][:, k]
    return f


def wall_fields(r):
    """The statistics of one wall_stats file (user_stats of tau, duprat_dpds
    and duprat_alpha): <|tau_w|>, <dp/ds> and <alpha>."""
    return {'tau': r['U'][:, 0], 'dpds': r['U'][:, 1], 'alpha': r['S'][:, 0]}


def combine(series, start, from_time=None, names=fields):
    """Read all windows and average them.

    Returns (geo, mean, windows, info): geo has the coordinates 'x', 'y'
    (nel, ly, lx) and 'lx'; mean is the length-weighted average of the
    windows used; windows is a list of (t_begin, t_end, fields); info lists
    the windows skipped and why. `names` turns a file into the dict of
    fields: `fields` for fluid_stats, `wall_fields` for wall_stats."""
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
            windows.append((a, b, names(r)))
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

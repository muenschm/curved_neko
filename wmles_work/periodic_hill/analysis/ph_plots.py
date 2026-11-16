#!/usr/bin/env python3
"""Figures of the periodic-hill statistics against the LES of Gloerfelt &
Cinnella.

Usage:
  python3 -I ph_plots.py --ref GLOERFELT_DIR
                         --run NAME DIR START [--run NAME DIR START ...]
                         [--from-time T] [--outdir DIR] [--tag TAG]
                         [--title TITLE]

GLOERFELT_DIR holds the database files Profiles_Re10595_grid*.dat and
Cf_Re10595_grid*.dat. DIR and START are as for ph_stats.py; NAME appears in
the legends. The runs get the colours and line styles in the order given
(at most 4). Written to --outdir, with _TAG appended to the names:

- profiles_u, profiles_v, profiles_upup, profiles_vpvp, profiles_upvp:
  profiles at x/h = 0.05, 0.5, 1, 2, ..., 8, drawn as x/h + scale * quantity
  over the hill (same scales as in the report);
- cf_bottom: C_f along the hill wall: the four reference grids, the
  experimental reattachment point, and for every run with wall statistics
  C_f of the Duprat model on the mean flow (solid) and sign(u_t) <|tau_w|>
  (dotted), see wall_summary.py;
- wall_alpha_dpds: <alpha> and <dp/ds> along the hill wall (runs with wall
  statistics only).
"""
import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ph_io import LX, LY, combine, hill_height, stats_files  # noqa: E402
from ph_io import wall_fields  # noqa: E402
from ph_stats import profiles  # noqa: E402
from wall_summary import cf_model, u_tangential, wall_nodes  # noqa: E402

XS = [0.05, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
# colour and line style per run, in this order (validated categorical
# palette; the line style is the second encoding)
STYLES = [('#2a78d6', '-'), ('#eb6834', '--'), ('#1baf7a', '-.'),
          ('#eda100', ':')]
REF_GRIDS = ['512x256x256', '256x128x128', '128x64x64', '64x33x32']
REF_GREYS = ['#000000', '#555555', '#999999', '#cccccc']
X_EXP = 4.21          # reattachment point of Rapp & Manhart
RE = 10595.0
JS = 3                # sampling row, as in the case file
# quantity: (column in the profile CSV of ph_stats, column in the reference
# file, sign, scale, axis label)
QUANT = {'u': (2, 1, 1.0, 0.8, r'x/h + 0.8 $\cdot\ \langle u\rangle/U_b$'),
         'v': (3, 2, 1.0, 5.0, r'x/h + 5 $\cdot\ \langle v\rangle/U_b$'),
         'upup': (4, 3, 1.0, 10.0,
                  r"x/h + 10 $\cdot\ \langle u'u'\rangle/U_b^2$"),
         'vpvp': (5, 4, 1.0, 20.0,
                  r"x/h + 20 $\cdot\ \langle v'v'\rangle/U_b^2$"),
         'upvp': (7, 5, -1.0, 20.0,
                  r"x/h + 20 $\cdot\ -\langle u'v'\rangle/U_b^2$")}


def read_zones(path):
    """The zones of a Tecplot-style ASCII file as a list of arrays."""
    zones, cur = [], None
    for line in open(path):
        s = line.strip()
        if s.startswith('ZONE'):
            cur = []
            zones.append(cur)
        elif cur is not None and s and (s[0].isdigit() or s[0] in '-.'):
            cur.append([float(v) for v in s.split()])
    return [np.array(z) for z in zones]


def hill_patch(ax):
    xs = np.linspace(0.0, LX, 901)
    hs = np.array([hill_height(x) for x in xs])
    ax.fill_between(xs, 0.0, hs, color='#d9d9d9', zorder=0)
    ax.plot(xs, hs, color='k', lw=0.8)
    ax.axhline(LY, color='k', lw=0.8)


def plot_profiles(q, ref, runs, title, path):
    col, rcol, sign, scale, label = QUANT[q]
    fig, ax = plt.subplots(figsize=(11, 4.2))
    hill_patch(ax)
    for x in XS:
        ax.axvline(x, color='#bbbbbb', lw=0.5, ls=':', zorder=0)
    for k, (x, z) in enumerate(zip(XS, ref)):
        ax.plot(x + scale*sign*z[:, rcol], z[:, 0], color='k', lw=1.0,
                label='Gloerfelt & Cinnella 512x256x256' if k == 0 else None)
    for (name, run), (c, ls) in zip(runs, STYLES):
        p = run['prof']
        for k, x in enumerate(XS):
            m = p[:, 0] == x
            ax.plot(x + scale*sign*p[m, col], p[m, 1], color=c, ls=ls,
                    lw=1.6, label=name if k == 0 else None)
    ax.set_xlim(-0.2, LX + 0.5)
    ax.set_ylim(0.0, LY + 0.05)
    ax.set_xlabel(label)
    ax.set_ylabel('y/h')
    ax.set_title(title, fontsize=10)
    ax.legend(loc='upper right', fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_cf(refdir, runs, title, path):
    fig, ax = plt.subplots(figsize=(10, 5))
    for g, grey in zip(REF_GRIDS, REF_GREYS):
        f = os.path.join(refdir, f'Cf_Re10595_grid{g}.dat')
        if os.path.exists(f):
            z = read_zones(f)[0]
            ax.plot(z[:, 0], z[:, 1], color=grey, lw=1.0,
                    label=f'Gloerfelt & Cinnella {g}')
    for (name, run), (c, ls) in zip(runs, STYLES):
        if run.get('cf_model') is None:
            continue
        ax.plot(run['xw'], run['cf_model'], color=c, ls=ls, lw=1.8,
                label=f'{name}, Duprat model on mean flow')
        ax.plot(run['xw'], run['cf_tau'], color=c, ls=':', lw=1.0,
                label=fr'{name}, sign$\cdot\langle|\tau_w|\rangle$')
    ax.axvline(X_EXP, color='m', ls=':', lw=1.0,
               label=f'exp. reattachment {X_EXP}')
    ax.axhline(0.0, color='k', lw=0.6)
    ax.set_xlim(0.0, LX)
    ax.set_xlabel('x/h')
    ax.set_ylabel(r'$C_f$')
    ax.set_title(title, fontsize=10)
    ax.legend(loc='upper left', fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_wall(runs, title, path):
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    for (name, run), (c, ls) in zip(runs, STYLES):
        if run.get('alpha') is None:
            continue
        a1.plot(run['xw'], run['alpha'], color=c, ls=ls, lw=1.6, label=name)
        a2.plot(run['xw'], run['dpds'], color=c, ls=ls, lw=1.6, label=name)
    a1.set_ylabel(r'$\langle\alpha\rangle$')
    a2.set_ylabel(r'$\langle\partial p/\partial s\rangle$')
    a2.axhline(0.0, color='k', lw=0.6)
    a2.set_xlabel('x/h (hill wall)')
    a2.set_xlim(0.0, LX)
    a1.set_title(title, fontsize=10)
    a1.legend(fontsize=8)
    a2.legend(fontsize=8, loc='lower left')
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def load_run(d, start, from_time):
    geo, mean, _, _ = combine(stats_files(d), start, from_time)
    run = {'prof': profiles(geo, mean, XS)}
    nodes = wall_nodes(geo, JS)
    run['xw'] = nodes['x']
    if glob.glob(os.path.join(d, 'wall_stats2d*.f[0-9]*')):
        _, wmean, _, _ = combine(stats_files(d, 'wall_stats2d'), start,
                                 from_time, names=wall_fields)
        ut = u_tangential(nodes, mean, JS)
        e, i = nodes['e'], nodes['i']
        run['cf_model'] = cf_model(nodes, ut, wmean, 1.0/RE)
        run['cf_tau'] = 2.0*np.sign(ut)*wmean['tau'][e, 0, i]
        run['alpha'] = wmean['alpha'][e, 0, i]
        run['dpds'] = wmean['dpds'][e, 0, i]
    return run


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--ref', required=True)
    ap.add_argument('--run', nargs=3, action='append', required=True,
                    metavar=('NAME', 'DIR', 'START'))
    ap.add_argument('--from-time', type=float, default=None)
    ap.add_argument('--outdir', default='.')
    ap.add_argument('--tag', default='')
    ap.add_argument('--title', default='')
    args = ap.parse_args()
    if len(args.run) > len(STYLES):
        raise SystemExit(f'at most {len(STYLES)} runs')
    runs = []
    for name, d, start in args.run:
        st = [float(v) for v in start.split(',')]
        runs.append((name, load_run(d, st, args.from_time)))
        print(f'read {name}')
    ref = read_zones(os.path.join(args.ref,
                                  'Profiles_Re10595_grid512x256x256.dat'))
    os.makedirs(args.outdir, exist_ok=True)
    tag = f'_{args.tag}' if args.tag else ''
    for q in QUANT:
        plot_profiles(q, ref, runs, args.title,
                      os.path.join(args.outdir, f'profiles_{q}{tag}.png'))
    plot_cf(args.ref, runs, args.title,
            os.path.join(args.outdir, f'cf_bottom{tag}.png'))
    plot_wall(runs, args.title,
              os.path.join(args.outdir, f'wall_alpha_dpds{tag}.png'))
    print(f'figures written to {args.outdir}')


if __name__ == '__main__':
    main()

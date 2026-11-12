# Periodic-hill post-processing

Scripts for the z-averaged statistics of the periodic-hill runs. They need
only Python 3 and numpy. Each script explains its inputs in the docstring at
the top; `--help` prints the options.

Run them with `python3 -I script.py ...` (or plain `python3`). The scripts
find each other through their own directory, so they can be called from
anywhere.

## Inputs

The case files write two kinds of 2D statistics (`"avg_direction": "z"`):

- `fluid_stats2d*.fNNNNN` from `fluid_stats` (`"set_of_stats": "basic"`):
  U = (<u>, <v>), P = <p>, T = <pp>, S = (<uu>, <vv>, <ww>, <uv>, <uw>, <vw>,
  <w>).
- `wall_stats2d*.fNNNNN` from `user_stats` of `tau`, `duprat_dpds` and
  `duprat_alpha`: U = (<tau>, <dp/ds>), S = (<alpha>). Duprat runs only.

Neko writes the average since the previous output and then resets the
statistics, so every file is one averaging window. The first file of a series
also holds the coordinates.

## Scripts

| script | what it does | example |
|---|---|---|
| `ph_stats.py` | mean velocity and Reynolds stresses along x = const; averages all windows (weighted by length), optionally from `--from-time` on; CSV output | `ph_stats.py RUN 150 --from-time 200 --out prof.csv` |
| `gpu_cpu_windows.py` | compares two runs window by window: peaks of u'u', -u'v', v'v' and min <u> at x/h 0.5...6, with standard errors from the scatter of the windows | `gpu_cpu_windows.py --run CPU dirA 150 --run GPU dirB 150` |
| `wall_summary.py` | separation and reattachment (sign change of the mean wall-parallel velocity at the sampling row), mean wall shear stress, alpha spikes at element boundaries | `wall_summary.py RUN1 RUN2` |
| `wall_dpds_check.py` | dp/ds used by the Duprat model against the gradient of the mean pressure, per wall region; also holds the file reader `read_fld` | `wall_dpds_check.py RUN1 RUN2` |

START is the time at which the statistics start (`start_time` of the stats
components). After a restart, give the restart time as well, e.g.
`150,483.3943`: a restarted run starts a new averaging window there.

## How the profiles are computed (`ph_stats.py`)

1. `read_fld` reads each file (Nek5000 format: header, element ids, then the
   blocks X, U, P, T, S; each scalar is a block over all elements).
2. `combine` turns the files into windows (t_begin, t_end, fields) and
   averages them.
3. `reynolds` forms the fluctuations, e.g. u'u' = <uu> - <u><u>.
4. `locate` finds, for every point of a profile, the element and the
   reference coordinates (r, s) with Newton's method on the element mapping;
   `interp` evaluates the GLL Lagrange interpolant there.

Points that lie just outside the discrete wall (where the hill meets the
flat floor inside an element) get NaN and are reported as "points outside".

## Checks against pySEMTools (`checks/`)

[pySEMTools](https://github.com/ExtremeFLOW/pySEMTools) (needs mpi4py) is
used as an independent check:

- `check_reader.py FILE ANALYSIS_DIR`: `read_fld` against `pynekread`. On
  12 Nov 2026 identical for all 11 fields of a `fluid_stats2d` file.
- `check_interp.py FILE PROFILES.csv OUTDIR`: the profiles of `ph_stats.py`
  against the pySEMTools probes. The probes need 3D elements, so the 2D
  data are extruded into one layer constant in z. On 12 Nov 2026: median
  difference 2e-16, 99th percentile 6e-15; only a single wall point that
  lies just outside the discrete wall differs (7e-4).

Note: until 12 Nov 2026 `read_fld` read multiple scalars in the wrong order
(element-major instead of one block per scalar). This affected only files
with more than one scalar, i.e. `fluid_stats2d`; the `wall_stats2d`
evaluations (one scalar) were correct.

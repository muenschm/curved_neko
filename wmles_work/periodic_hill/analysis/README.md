# Periodic-hill post-processing

Scripts for the z-averaged statistics of the periodic-hill runs. They need
only Python 3 and numpy (`ph_plots.py` also matplotlib). Each script explains its inputs in the docstring at
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
| `ph_stats.py` | mean velocity and Reynolds stresses along x = const; averages all windows (weighted by length), optionally from `--from-time` on; CSV output | `ph_stats.py RUN 150 --x 0.5,1,2,3,4,6 --out prof.csv` |
| `wall_summary.py` | separation and reattachment (sign of the mean wall-parallel velocity at the sampling point, and of C_f reconstructed with the Duprat model), C_f extrema, mean wall shear stress, reattachment per window with standard error | `wall_summary.py RUN 150 --out wall.csv` |
| `wall_dpds_check.py` | dp/ds used by the Duprat model against the gradient of the mean pressure, per wall region | `wall_dpds_check.py RUN 150` |
| `gpu_cpu_windows.py` | compares two runs window by window: peaks of u'u', -u'v', v'v' and min <u> at x/h 0.5...6, with standard errors from the scatter of the windows | `gpu_cpu_windows.py --run CPU dirA 150 --run GPU dirB 150` |
| `ph_forces.py` | time-averaged forces on both walls, forcing and momentum balance, from the CSV files of `../forces_from_log.py` | `ph_forces.py FORCEDIR --from-time 150` |
| `ph_plots.py` | figures as in the report: profiles of <u>, <v>, u'u', v'v', -u'v' over the hill, C_f along the hill wall with the four reference grids, <alpha> and <dp/ds> along the wall; up to 4 runs | `ph_plots.py --ref GLOERFELT_DIR --run "new, Fritz" dirA 150 --run "old, Fritz" dirB 150,483.3943 --outdir figs --tag v2` |

Modules used by the scripts:

- `ph_io.py`: hill geometry, the file reader `read_fld`, the grouping of the
  files into series and the averaging of the windows (`combine`).
- `duprat_model.py`: the Duprat wall model in numpy (same equations as the
  Neko kernel, including the reversed wall shear stress); run it directly to
  compare it with the reference values of the unit tests.

START is the time at which the statistics start (`start_time` of the stats
components). After a restart, give the restart time as well, e.g.
`150,489.53`: a restarted run starts a new averaging window there. The
restart time is the time of the checkpoint, i.e. the last time of the
previous job.

Put the `fluid_stats2d*` and `wall_stats2d*` files of a run (both jobs, with
the `.nek5000` files) into one directory.

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

- `check_reader.py FILE ANALYSIS_DIR`: `read_fld` (ph_io.py) against `pynekread`. On
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

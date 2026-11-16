# Periodic hill, Duprat, version bcb0e5a (paper version)

Production runs of `feature/duprat-local-dpds` at `bcb0e5a`: no Stokes limit
(`stokes_limit` false), reversed wall shear stress, gradient averaged across
element boundaries and filtered as a vector (`filter_time` 1, from t = 50).
N = 8, statistics from t = 150, windows of 50 time units, end t = 600.

- Fritz (CPU, 72 ranks): two jobs, restart at t = 489.5319 (series
  `*2d00` and `*2d60`).
- Alex (GPU, CUDA, 1 rank): one job.

Files:
- `summary.txt`: output of all scripts (below).
- `wall_fritz.csv`, `wall_alex.csv`: along the hill wall, x, y, <u_t>, h,
  <dp/ds>, <alpha>, C_f (Duprat model from the mean flow), C_f (from
  <|tau_w|>, signed with <u_t>).
- `profiles_fritz.csv`, `profiles_alex.csv`: profiles at x/h = 0.5, 1, 2, 3,
  4, 6.
- `gpu_cpu_windows.txt`: Fritz against Alex, window by window.

Commands (run directory with all stats files of a run; force CSVs from
`../forces_from_log.py` on the logs):

```bash
A=../analysis
python3 -I $A/wall_summary.py    fritz 150,489.5319246170351 --out wall_fritz.csv
python3 -I $A/wall_summary.py    alex  150 --out wall_alex.csv
python3 -I $A/wall_dpds_check.py fritz 150,489.5319246170351
python3 -I $A/wall_dpds_check.py alex  150
python3 -I $A/ph_stats.py        fritz 150,489.5319246170351 --x 0.5,1,2,3,4,6 --out profiles_fritz.csv
python3 -I $A/ph_stats.py        alex  150 --x 0.5,1,2,3,4,6 --out profiles_alex.csv
python3 -I $A/gpu_cpu_windows.py --run Fritz fritz 150,489.5319246170351 --run Alex alex 150
python3 ../forces_from_log.py log_1Node_1st_run log_1Node_2nd_run   # in forces_fritz/
python3 -I $A/ph_forces.py       forces_fritz
```

## Comparison with the previous Duprat version and the reference

Old Duprat (`feature/duprat-wall-model`) and reference values from the
report (Tables hilldupsep, hilldupforces, hilldupprof). Fritz / Alex.

| | old Duprat | new Duprat (bcb0e5a) | reference |
|---|---|---|---|
| separation x/h | 0.45 / 0.45 | 0.44 / 0.44 | -0.12 (LES), 0.19 (Breuer) |
| reattachment x/h, mean field | 3.33 / 3.21 | 3.46 / 3.64 | 4.31 (LES), 4.21 (exp.) |
| reattachment per window | 3.37 +- 0.10 / 3.21 +- 0.08 | 3.43 +- 0.11 / 3.64 +- 0.04 | |
| hill pressure drag Fx_p | 0.0622 / 0.0595 | 0.0568 / 0.0559 | |
| hill skin friction Fx_v | 0.00285 / 0.00286 | 0.00253 / 0.00234 | |
| top skin friction Fx_v | 0.00854 / 0.00859 | 0.00860 / 0.00874 | |
| forcing <f> | 0.01291 / 0.01243 | 0.01190 / 0.01174 | |
| C_f min (bubble), model from mean flow | -0.0015 / -0.0013 (*) | -0.0045 / -0.0045 | -0.0038...-0.0044 |
| C_f max (windward) | 0.031 / 0.031 (*) | 0.038 / 0.036 | 0.030...0.033 |

(*) reconstructed in the report with the previous model.

Model dp/ds against the gradient of the mean pressure (ratio of the region
means): old about 0.55 (crest), 1/13 (bubble), 1/8 (recovery), 2/3
(windward); new 0.85, 0.74 / 0.81, 0.24, 1.00.

Profiles (min <u>, max u'u', max -u'v'), Fritz / Alex:

| x/h | min <u> old | new | ref | u'u' old | new | ref | -u'v' old | new | ref |
|---|---|---|---|---|---|---|---|---|---|
| 0.5 | -0.112 / -0.111 | -0.109 / -0.111 | -0.080 | 0.129 / 0.127 | 0.117 / 0.113 | 0.112 | 0.046 / 0.047 | 0.044 / 0.042 | 0.042 |
| 1 | -0.149 / -0.137 | -0.124 / -0.126 | -0.135 | 0.110 / 0.108 | 0.100 / 0.097 | 0.093 | 0.051 / 0.052 | 0.049 / 0.046 | 0.042 |
| 2 | -0.237 / -0.197 | -0.209 / -0.211 | -0.253 | 0.104 / 0.103 | 0.089 / 0.086 | 0.087 | 0.046 / 0.044 | 0.040 / 0.037 | 0.036 |
| 3 | -0.115 / -0.084 | -0.095 / -0.119 | -0.140 | 0.100 / 0.099 | 0.091 / 0.085 | 0.091 | 0.041 / 0.039 | 0.034 / 0.033 | 0.035 |
| 4 | slip | +0.005 / -0.013 | -0.005 | 0.086 / 0.087 | 0.079 / 0.078 | 0.079 | 0.034 / 0.033 | 0.028 / 0.029 | 0.030 |
| 6 | slip | +0.084 / +0.084 | 0.000 | 0.064 / 0.064 | 0.062 / 0.059 | 0.055 | 0.023 / 0.022 | 0.021 / 0.020 | 0.020 |

Notes:
- The mean of <|tau_w|> signed with <u_t> is 16 % above the force_torque
  friction of the hill wall (old run: 0.4 %). The reversed branch acts
  instantaneously: tau_w points backwards at times while the mean flow at the
  sampling point goes forwards. Use force_torque for the friction.
- In the mean flow, the model never gives tau_w against <u_t>, so the zero
  crossings of <u_t> and of the reconstructed C_f agree.
- Fritz against Alex: 17 of 18 turbulence peaks higher on the CPU, +2 to
  +10 %, all below 2 sigma; the same tendency as Reichardt (16 of 18).

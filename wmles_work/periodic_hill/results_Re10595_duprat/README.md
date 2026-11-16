# Periodic hill, Duprat, previous version (feature/duprat-wall-model)

Re-evaluation with the scripts in `../analysis` (16 Nov 2026) of the run
described in the report: Stokes limit on, projection on the instantaneous
flow direction before the time filter, no reversed wall shear stress.

- Fritz (CPU, 72 ranks): job 1 to t = 493.75 with statistics up to t = 450
  (the window 450-493.75 was not written), job 2 restarted at t = 483.3943.
  Windows 150-450 and 483.4-600; START = 150,483.3943.
- Alex (GPU, CUDA, 1 rank): one job, windows 150-600; START = 150.

The scripts reproduce the numbers of the report for both runs: separation
0.45 / 0.45, reattachment 3.33 / 3.21 (Alex per window 3.21 +- 0.08), C_f
min -0.0015 / -0.0013 and max 0.031, the dp/ds ratios, the profile table,
and the mean C_f from <|tau_w|> 0.00284 (force_torque 0.00285 / 0.00286).
For Fritz, the reattachment per window is 3.40 +- 0.11 here against
3.37 +- 0.10 in the report; here windows shorter than 40 (483.4-500) are
left out.

Files: `summary_fritz.txt`, `summary_alex.txt`, `wall_*.csv`,
`profiles_*.csv`, `gpu_cpu_windows.txt` (Fritz against Alex).

```bash
A=../analysis
python3 -I $A/wall_summary.py    old_fritz 150,483.3943 --out wall_fritz.csv
python3 -I $A/wall_dpds_check.py old_fritz 150,483.3943
python3 -I $A/ph_stats.py        old_fritz 150,483.3943 --x 0.5,1,2,3,4,6 --out profiles_fritz.csv
python3 -I $A/wall_summary.py    old_alex 150 --out wall_alex.csv
python3 -I $A/wall_dpds_check.py old_alex 150
python3 -I $A/ph_stats.py        old_alex 150 --x 0.5,1,2,3,4,6 --out profiles_alex.csv
python3 -I $A/gpu_cpu_windows.py --run old_Fritz old_fritz 150,483.3943 --run old_Alex old_alex 150
```

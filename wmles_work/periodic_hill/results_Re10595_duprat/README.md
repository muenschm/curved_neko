# Periodic hill, Duprat, previous version (feature/duprat-wall-model)

Re-evaluation with the scripts in `../analysis` (16 Nov 2026) of the run
described in the report: Stokes limit on, projection on the instantaneous
flow direction before the time filter, no reversed wall shear stress.

- Fritz (CPU, 72 ranks): job 1 to t = 493.75 with statistics up to t = 450
  (the window 450-493.75 was not written), job 2 restarted at t = 483.3943.
  Windows 150-450 and 483.4-600; START = 150,483.3943.
- Alex: to follow.

The scripts reproduce the numbers of the report for this run: separation
0.45, reattachment 3.33, C_f min -0.0015 / max 0.031, the dp/ds ratios, the
profile table, and the mean C_f from <|tau_w|> 0.00284 (force_torque
0.00285).

```bash
A=../analysis
python3 -I $A/wall_summary.py    old_fritz 150,483.3943 --out wall_fritz.csv
python3 -I $A/wall_dpds_check.py old_fritz 150,483.3943
python3 -I $A/ph_stats.py        old_fritz 150,483.3943 --x 0.5,1,2,3,4,6 --out profiles_fritz.csv
```

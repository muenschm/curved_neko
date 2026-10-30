#!/usr/bin/env python3
"""Extract the force_torque and flow-rate output from Neko log files into CSV.

Usage:  python3 forces_from_log.py neko.log [more.log ...]

Writes one file per zone, forces_zone<id>.csv, with the columns
step, time, Fx, Fx_p, Fx_v, Fy, Fy_p, Fy_v, Fz, Fz_p, Fz_v
(total, pressure and viscous part of each force component, as scaled by the
"scale" of the force_torque block), and forcing.csv with the columns
step, time, f, where f is the body force per unit volume applied by
flow_rate_force (the "Scale" column of the "Flow rate" log lines, i.e. the
mean pressure gradient -dP/dx that drives the flow). Several logs, e.g. of
restarted runs, are concatenated and sorted by time.

For the periodic hill case, scale = 2 / (rho U_b^2 L_x L_z), so Fx_p and Fx_v
are the pressure-drag and mean skin-friction coefficients of the wall. With
the flow rate held constant, the time average of f * V equals the time
average of the total x force on both walls, (Fx_bottom + Fx_top) / scale.
"""
import re
import sys

HEADER = re.compile(r"Force and torque on zone\s+(\d+)")
LINE = re.compile(r"^\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+),\s*force([xyz])\s*$")
FLOW = re.compile(r"^\s*(\d+)\s+\|\s+[xyz]\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s*$")


def parse(paths):
    forces = {}  # zone -> {(step, time): {component: (total, p, v)}}
    forcing = {}  # (step, time) -> f
    for path in paths:
        zone = None
        with open(path, errors="ignore") as f:
            for line in f:
                m = HEADER.search(line)
                if m:
                    zone = int(m.group(1))
                    continue
                m = LINE.match(line)
                if m and zone is not None:
                    step, t = int(m.group(1)), float(m.group(2))
                    vals = tuple(float(m.group(k)) for k in (3, 4, 5))
                    forces.setdefault(zone, {}).setdefault((step, t), {})[m.group(6)] = vals
                    continue
                m = FLOW.match(line)
                if m:
                    forcing[(int(m.group(1)), float(m.group(2)))] = float(m.group(3))
    return forces, forcing


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    forces, forcing = parse(sys.argv[1:])
    if not forces:
        sys.exit("No force_torque output found.")
    for zone, rows in sorted(forces.items()):
        name = f"forces_zone{zone}.csv"
        n = 0
        with open(name, "w") as out:
            out.write("step,time,Fx,Fx_p,Fx_v,Fy,Fy_p,Fy_v,Fz,Fz_p,Fz_v\n")
            for (step, t), comp in sorted(rows.items(), key=lambda r: r[0][1]):
                if len(comp) < 3:
                    continue
                vals = [v for c in "xyz" for v in comp[c]]
                out.write(f"{step},{t:.10e}," + ",".join(f"{v:.10e}" for v in vals) + "\n")
                n += 1
        print(f"zone {zone}: {n} samples -> {name}")
    if forcing:
        with open("forcing.csv", "w") as out:
            out.write("step,time,f\n")
            for (step, t), f in sorted(forcing.items(), key=lambda r: r[0][1]):
                out.write(f"{step},{t:.10e},{f:.10e}\n")
        print(f"flow-rate forcing: {len(forcing)} samples -> forcing.csv")


if __name__ == "__main__":
    main()

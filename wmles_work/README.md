# WMLES work files

Working files of the wall-model work (Reichardt and Duprat) for Neko: case
files, evaluation scripts, results, report sources and pull request texts.

The cloud sessions in which this work is done run in short-lived containers.
To keep the files, `tools/autosave.sh` commits and pushes this directory to the
branch `claude/stoic-euler-8zyrhx` after every turn. A Claude Code Stop hook in
`~/.claude/settings.json` calls it.

- **Branch:** https://github.com/muenschm/curved_neko/tree/claude/stoic-euler-8zyrhx/wmles_work
- **Large data:** simulation output (`*.f00000`, checkpoints, meshes, logs,
  archives) is not committed; see `.gitignore`. The script also skips any file
  above 20 MB.

## Contents

- `pr_texts/`: descriptions and comments for the upstream pull requests
  - `upstream_pr_reichardt.md`: [ExtremeFLOW/neko#2884](https://github.com/ExtremeFLOW/neko/pull/2884)
  - `upstream_pr_reichardt_comment.md`: comment on that pull request after the
    develop merge and the removal of the GPU warning
  - `upstream_pr_duprat.md`: [ExtremeFLOW/neko#2888](https://github.com/ExtremeFLOW/neko/pull/2888)
- `report/`: the WMLES report, `wmles_report.tex` with `figures/`, and the
  PDF built from it (`pdflatex` twice). It includes the periodic-hill results
  for Reichardt and Duprat and the removal of the GPU warnings.
- `periodic_hill/`: the periodic-hill case and its evaluation.
  - `forces_from_log.py`: extracts the `force_torque` and flow-rate output
    from Neko logs into CSV.
  - `analysis/gpu_cpu_windows.py`: compares two runs window by window. It
    imports `ph_stats.py`.
- `tools/autosave.sh`: the autosave script.

### Lost and still to be restored

The container reset on 30 Oct 2026 deleted these files before they were
saved here:
- the case files (`periodic_hill.case`, the Duprat case, `periodic_hill.f90`,
  `CHANGES.md`, `gpu_cpu_check/`);
- the evaluation scripts `analysis/ph_stats.py` and `analysis/ph_forces.py`,
  with the results folders.

They are rewritten here when they are next needed.

## Restoring the hook in a new session

If a new container lacks `~/.claude/settings.json`, add this Stop hook:

```json
{
  "hooks": {
    "Stop": [
      {
        "hooks": [
          {
            "type": "command",
            "command": "/home/user/curved_neko/wmles_work/tools/autosave.sh",
            "timeout": 60
          }
        ]
      }
    ]
  }
}
```

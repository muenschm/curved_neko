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
- `tools/autosave.sh`: the autosave script.

The periodic-hill case, the evaluation scripts and the report are added here
once they are restored.

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

#!/usr/bin/env bash
# Commit and push the WMLES work files (wmles_work/) on the session branch.
# Runs from a Claude Code Stop hook after every turn, so that nothing is lost
# when the cloud container is reset. It never blocks the session: on any
# problem it reports and exits 0.
set -u
REPO=/home/user/curved_neko
BRANCH=claude/stoic-euler-8zyrhx
MAX_BYTES=$((20 * 1024 * 1024))

cd "$REPO" 2>/dev/null || exit 0
[ "$(git rev-parse --abbrev-ref HEAD 2>/dev/null)" = "$BRANCH" ] || exit 0

git add -A -- wmles_work >/dev/null 2>&1 || exit 0
git diff --cached --quiet -- wmles_work && exit 0

# Leave out files that are too large for the repository
git diff --cached --name-only -z --diff-filter=AM -- wmles_work |
  while IFS= read -r -d '' f; do
    if [ -f "$f" ] && [ "$(stat -c %s "$f")" -gt "$MAX_BYTES" ]; then
      git reset -q -- "$f"
    fi
  done
git diff --cached --quiet -- wmles_work && exit 0

git commit -q -m "Autosave WMLES work files ($(date -u +%Y-%m-%dT%H:%MZ))

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01Qr9jBZiogCtaH7KGcZPRSx" || exit 0

for wait in 0 2 4 8; do
  sleep "$wait"
  if git push -q -u origin "$BRANCH" >/dev/null 2>&1; then
    exit 0
  fi
done
echo '{"systemMessage": "Autosave: the push of wmles_work/ failed; the commit is kept locally."}'
exit 0

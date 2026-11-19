#!/usr/bin/env bash
# Shows git state and the size of files git would track.
cd "$HOME/ngs-platform"
echo "== git version: $(git --version)"
echo "== user: $(git config --global user.name) / email set: $([ -n "$(git config --global user.email)" ] && echo yes || echo no)"
if [ -d .git ]; then
  echo "== repo exists, branch: $(git branch --show-current)"
  git log --oneline | head -5
  git status --short | head -30
else
  echo "== no repo yet"
fi
echo "== files that would be tracked (size) =="
git ls-files --others --cached --exclude-standard 2>/dev/null | xargs -r du -ch 2>/dev/null | tail -1 || \
  find . -path ./.git -prune -o -type f -print | xargs du -ch | tail -1
echo "== nextflow running? =="
pgrep -fa "nextflow" | grep -v pgrep | head -3 || echo "not running"
ls -d pipeline/work pipeline/.nextflow 2>/dev/null

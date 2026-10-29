#!/usr/bin/env bash
# Explicit manual publication step. Usage: bash scripts/publish_github.sh OWNER/REPO
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
repo="${1:?Specify OWNER/REPO}"
[[ "$repo" =~ ^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$ ]] || { echo 'Use OWNER/REPO'; exit 1; }
command -v gh >/dev/null
gh auth status
if [[ -e .git ]]; then
  echo 'This directory already has Git metadata. Use the git commands in README to avoid changing an existing remote.'
  exit 1
fi
git init -b main
git add .
git commit -m 'Add reproducible RRBS analysis workflow'
gh repo create "$repo" --private --source=. --remote=origin --push

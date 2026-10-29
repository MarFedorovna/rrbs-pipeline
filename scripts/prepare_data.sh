#!/usr/bin/env bash
# Download reads and reference annotations, then build the bisulfite index.
# Run inside the activated environment. Arguments are passed to Snakemake.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
bash scripts/run.sh prepare_data "$@"

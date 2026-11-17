#!/usr/bin/env bash
# Runs the pipeline on the nf-core test reads.
# Usage: bash pipeline/run_test.sh [test|test_grch38] [-resume]
#   test        : 40 kb chr22 fragment reference (1 minute, checks the pipeline runs)
#   test_grch38 : full GRCh38 reference, calling on chr22 (10-20 minutes, finds real variants)
set -euo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ngs
export PATH="$HOME/.local/bin:$PATH"
cd "$(dirname "$0")"
MODE=${1:-test_grch38}
shift || true
if [ "$MODE" = "test" ]; then
  PROFILE="test,conda"
else
  PROFILE="grch38,test_grch38,conda"
fi
nextflow run main.nf -profile "$PROFILE" -ansi-log false "$@"

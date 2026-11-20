#!/usr/bin/env bash
# Creates the 'wesapp' Python env, installs the interpretation package and downloads
# ClinVar's variant_summary (used to build the PS1/PM5 protein-change index).
# Usage: bash scripts/setup/06_setup_interpretation.sh
set -euo pipefail
REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
source "$HOME/miniforge3/etc/profile.d/conda.sh"

if conda env list | grep -qE '^wesapp\s'; then
  mamba env update -n wesapp -f "$REPO_DIR/envs/wesapp.yml"
else
  mamba env create -f "$REPO_DIR/envs/wesapp.yml"
fi

D="$HOME/ngs-data/annotation/clinvar"
mkdir -p "$D"
echo ">> Downloading ClinVar variant_summary..."
wget -c -q -O "$D/variant_summary.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
ls -lh "$D/variant_summary.txt.gz"
echo ">> Done"

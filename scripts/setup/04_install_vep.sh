#!/usr/bin/env bash
# Creates the 'vep' conda environment and downloads the VEP plugins we use.
# Usage:  bash scripts/setup/04_install_vep.sh
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
RELEASE=116
PLUGINS_DIR="$HOME/ngs-data/annotation/vep_plugins"

source "$HOME/miniforge3/etc/profile.d/conda.sh"

if conda env list | grep -qE '^vep\s'; then
  echo ">> 'vep' environment exists, updating..."
  mamba env update -n vep -f "$REPO_DIR/envs/vep.yml" --prune
else
  echo ">> Creating 'vep' environment (5-10 minutes)..."
  mamba env create -f "$REPO_DIR/envs/vep.yml"
fi

echo ">> Downloading VEP plugins (release/$RELEASE)..."
mkdir -p "$PLUGINS_DIR"
for p in AlphaMissense REVEL; do
  curl -fsSL -o "$PLUGINS_DIR/$p.pm" \
    "https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/$RELEASE/$p.pm"
done
ls -l "$PLUGINS_DIR"

conda activate vep
echo ">> VEP version:"
vep --help 2>&1 | grep -E 'ensembl-vep|ensembl ' | head -3

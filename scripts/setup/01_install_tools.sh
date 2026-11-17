#!/usr/bin/env bash
# Installs Miniforge (conda/mamba), creates the 'ngs' tool environment and installs Nextflow.
# No sudo needed: everything goes into your home folder.
# Usage:  bash scripts/setup/01_install_tools.sh
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"

# ---------- 1. Miniforge ----------
if [ ! -d "$HOME/miniforge3" ]; then
  echo ">> Installing Miniforge..."
  wget -q --show-progress -O /tmp/miniforge.sh \
    https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
  rm /tmp/miniforge.sh
  "$HOME/miniforge3/bin/conda" init bash
else
  echo ">> Miniforge already installed, skipping."
fi

# Make conda/mamba usable inside this script
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda config --set channel_priority strict

# ---------- 2. 'ngs' tool environment ----------
if conda env list | grep -qE '^ngs\s'; then
  echo ">> Updating existing 'ngs' environment..."
  mamba env update -n ngs -f "$REPO_DIR/envs/ngs-tools.yml" --prune
else
  echo ">> Creating 'ngs' environment (takes 5-15 minutes)..."
  mamba env create -f "$REPO_DIR/envs/ngs-tools.yml"
fi

# ---------- 3. Nextflow ----------
mkdir -p "$HOME/.local/bin"
if [ ! -x "$HOME/.local/bin/nextflow" ]; then
  echo ">> Installing Nextflow..."
  (cd "$HOME/.local/bin" && curl -s https://get.nextflow.io | bash)
  chmod +x "$HOME/.local/bin/nextflow"
fi
export PATH="$HOME/.local/bin:$PATH"

# ---------- 4. Version check ----------
conda activate ngs
echo
echo "================ Installed versions ================"
echo "nextflow : $(nextflow -v 2>/dev/null | tr -s ' ')"
echo "bwa      : $(bwa 2>&1 | grep -i '^version' | cut -d' ' -f2)"
echo "samtools : $(samtools --version | head -1)"
echo "bcftools : $(bcftools --version | head -1)"
echo "gatk     : $(gatk --version 2>/dev/null | grep -i 'GATK' | head -1)"
echo "fastp    : $(fastp --version 2>&1)"
echo "mosdepth : $(mosdepth --version)"
echo "===================================================="
echo "Done. Close and reopen the terminal, then run:  conda activate ngs"

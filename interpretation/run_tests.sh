#!/usr/bin/env bash
# Installs the package in editable mode (first time) and runs the unit tests.
set -euo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate wesapp
cd "$(dirname "$0")"
pip show wesinterp >/dev/null 2>&1 || pip install -q -e .
python -m pytest -q "$@"

#!/usr/bin/env bash
# Adds the PDF reporting libraries (Jinja2, WeasyPrint, fonts) to the 'wesapp' environment.
set -euo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"
mamba install -y -q -n wesapp -c conda-forge jinja2 weasyprint fonts-conda-ecosystem
conda activate wesapp
python -c "import weasyprint, jinja2; print('weasyprint', weasyprint.__version__, '| jinja2', jinja2.__version__)"
echo "fonts available:"; fc-list : family | sort -u | head -8

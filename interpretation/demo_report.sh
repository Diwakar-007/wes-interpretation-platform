#!/usr/bin/env bash
# Demo: report for the simulated cystic fibrosis patient (known-variant VCF), with QC from the chr21 test run.
set -euo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate wesapp
OUT="$HOME/ngs-data/runs/demo_report"
mkdir -p "$OUT"
H="$HOME/ngs-data/annotation/hpo"
TERMS=$(for q in "bronchiectasis" "exocrine pancreatic insufficiency" "elevated sweat chloride"; do
  wesinterp hpo-search "$q" --hpo-dir "$H" --limit 1 | cut -f1; done | paste -sd,)
bash "$(dirname "$0")/run_classify.sh" /tmp/known_check.vep.vcf "$OUT/DEMO_CF" "$TERMS" 2>/dev/null
wesinterp report --json "$OUT/DEMO_CF.acmg.json" \
  --out-html "$OUT/DEMO_CF.report.html" --out-pdf "$OUT/DEMO_CF.report.pdf" \
  --sample-id TEST01 --sex "Not provided" \
  --indication "Recurrent chest infections and poor weight gain since infancy; elevated sweat chloride." \
  --sample-dir "$HOME/ngs-data/runs/test_grch38/TEST01"
ls -lh "$OUT"

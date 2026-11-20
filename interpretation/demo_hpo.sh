#!/usr/bin/env bash
# Demo: search HPO terms, then classify + rank the known-variant test VCF for a cystic-fibrosis-like patient.
set -euo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate wesapp
H="$HOME/ngs-data/annotation/hpo"
for q in "bronchiectasis" "pancreatic insufficiency" "elevated sweat chloride"; do
  echo "== hpo-search '$q'"; wesinterp hpo-search "$q" --hpo-dir "$H" --limit 3
done
TERMS=$(for q in "bronchiectasis" "exocrine pancreatic insufficiency" "elevated sweat chloride"; do
  wesinterp hpo-search "$q" --hpo-dir "$H" --limit 1 | cut -f1; done | paste -sd,)
echo; echo "== patient terms: $TERMS"; echo
time bash "$(dirname "$0")/run_classify.sh" /tmp/known_check.vep.vcf /tmp/known_check "$TERMS"
echo
cut -f6,10,12,15,18-21 /tmp/known_check.acmg.tsv | column -t -s $'\t'
echo
column -t -s $'\t' /tmp/known_check.candidates.tsv

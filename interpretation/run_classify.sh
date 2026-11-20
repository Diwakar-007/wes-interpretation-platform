#!/usr/bin/env bash
# Builds the ClinVar protein index (first time only) and classifies a VEP-annotated VCF.
# Usage: bash interpretation/run_classify.sh <sample.vep.vcf.gz> <out_prefix>
set -euo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate wesapp
A="$HOME/ngs-data/annotation"
IDX="$A/clinvar/clinvar_protein_index.tsv.gz"
VCF=${1:?VEP VCF required}
OUT=${2:?output prefix required}

if [ ! -s "$IDX" ]; then
  wesinterp build-clinvar-index --variant-summary "$A/clinvar/variant_summary.txt.gz" --out "$IDX"
fi

wesinterp classify --vcf "$VCF" \
  --loeuf "$A/gnomad/gnomad.v4.1.constraint_metrics.tsv" \
  --clinvar-index "$IDX" \
  --out "$OUT.acmg.tsv" --json "$OUT.acmg.json"

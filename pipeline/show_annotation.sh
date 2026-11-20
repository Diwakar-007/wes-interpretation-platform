#!/usr/bin/env bash
# Prints a readable table of VEP annotations (one PICKed consequence per variant per gene).
# Usage: bash pipeline/show_annotation.sh <sample.vep.vcf.gz> [n_rows]
set -euo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate vep
VCF=${1:-$HOME/ngs-data/runs/test_grch38/TEST01/vcf/TEST01.vep.vcf.gz}
N=${2:-25}

echo "== CSQ fields added by VEP =="
bcftools +split-vep -l "$VCF" | cut -f2 | tr '\n' ' ' | fold -w 150; echo; echo

echo "== consequence summary (PASS, picked) =="
bcftools +split-vep -f '%Consequence\n' -d -i 'FILTER="PASS" && PICK="1"' "$VCF" | sort | uniq -c | sort -nr
echo

echo "== variants (PASS, picked) =="
{
  echo -e "POS\tREF>ALT\tGT\tGENE\tCONSEQUENCE\tHGVSc\tHGVSp\tMANE\trsID\tgnomADe_AF\tgnomADg_AF\tClinVar\tREVEL\tAlphaMissense"
  bcftools +split-vep -d -i 'FILTER="PASS" && PICK="1"' \
    -f '%POS\t%REF>%ALT\t[%GT]\t%SYMBOL\t%Consequence\t%HGVSc\t%HGVSp\t%MANE_SELECT\t%Existing_variation\t%gnomADe_AF\t%gnomADg_AF\t%ClinVar_CLNSIG\t%REVEL\t%am_class\n' "$VCF" \
  | awk 'BEGIN{FS=OFS="\t"} {sub(/^[^:]*:/,"",$6); sub(/^[^:]*:/,"",$7); print}' \
  | head -n "$N"
} | column -t -s $'\t'

#!/usr/bin/env bash
# Quick summary of a finished run.  Usage: bash pipeline/summarise_run.sh <outdir> <sample>
set -euo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ngs
OUT=${1:-$HOME/ngs-data/runs/test}
S=${2:-TEST01}
VCF="$OUT/$S/vcf/$S.filtered.vcf.gz"

echo "== files =="
find "$OUT/$S" -type f | sed "s#$OUT/##" | sort
echo
echo "== mapping (samtools flagstat) =="
samtools flagstat "$OUT/$S/bam/$S.recal.bam" | grep -E "in total|primary mapped|properly paired|duplicates" | head -5
echo
echo "== duplication (MarkDuplicates) =="
grep -A1 "^LIBRARY" "$OUT/$S/qc/markduplicates/$S.md_metrics.txt" | cut -f9 | paste - -
echo
echo "== coverage over targets (mosdepth) =="
grep -E "total_region" "$OUT/$S/qc/mosdepth/$S.mosdepth.summary.txt"
echo
echo "== variants =="
echo "total records : $(bcftools view -H "$VCF" | wc -l)"
echo "PASS          : $(bcftools view -H -f PASS "$VCF" | wc -l)"
echo "PASS SNPs     : $(bcftools view -H -f PASS -v snps "$VCF" | wc -l)"
echo "PASS indels   : $(bcftools view -H -f PASS -v indels "$VCF" | wc -l)"
echo "FILTER values :"; bcftools query -f '%FILTER\n' "$VCF" | sort | uniq -c
echo
echo "== first 5 PASS variants =="
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT\t%DP]\n' -i 'FILTER="PASS"' "$VCF" | head -5

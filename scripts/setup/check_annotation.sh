#!/usr/bin/env bash
# Verifies the annotation resources: size, index, and a sample lookup in each.
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate vep
A="$HOME/ngs-data/annotation"
echo "== sizes / done markers =="
for d in vep_cache clinvar gnomad alphamissense revel vep_plugins; do
  printf "%-14s %8s  %s\n" "$d" "$(du -sh "$A/$d" 2>/dev/null | cut -f1)" "$([ -f "$A/$d/.done" ] && echo done || echo '-')"
done
echo
echo "== ClinVar: $(cat "$A/clinvar/VERSION.txt" 2>/dev/null), contigs: $(tabix -l "$A/clinvar/clinvar.chr.vcf.gz" | head -3 | tr '\n' ' ')"
tabix "$A/clinvar/clinvar.chr.vcf.gz" chr17:43044295-43125483 | head -2 | cut -f1-5,8 | cut -c1-200
echo
echo "== AlphaMissense header + lookup (chr21) =="
zcat "$A/alphamissense/AlphaMissense_hg38.tsv.gz" 2>/dev/null | sed -n '4p'
tabix "$A/alphamissense/AlphaMissense_hg38.tsv.gz" chr21:14965206-14975206 | head -2
echo
echo "== REVEL header + lookup (chr21) =="
zcat "$A/revel/revel_grch38_chr.tsv.gz" 2>/dev/null | head -1
tabix "$A/revel/revel_grch38_chr.tsv.gz" chr21:15819233-15919233 | head -2
echo
echo "== gnomAD constraint columns (first 12) =="
head -1 "$A/gnomad/gnomad.v4.1.constraint_metrics.tsv" | tr '\t' '\n' | head -12 | tr '\n' ' '; echo
echo "rows: $(wc -l < "$A/gnomad/gnomad.v4.1.constraint_metrics.tsv")"

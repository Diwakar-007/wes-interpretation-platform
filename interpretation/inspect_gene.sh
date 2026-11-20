#!/usr/bin/env bash
# Shows a gene's HPO inheritance terms per disease, and its number of phenotype annotations.
# Usage: bash interpretation/inspect_gene.sh CFTR
H="$HOME/ngs-data/annotation/hpo"
G=${1:-CFTR}
echo "== $G inheritance terms by disease =="
awk -F'\t' -v g="$G" '$2==g && $4 ~ /inheritance/ {print $6"\t"$4}' "$H/genes_to_phenotype.txt" | sort -u
echo "== diseases (genes_to_disease) =="
awk -F'\t' -v g="$G" '$2==g {print $3"\t"$4}' "$H/genes_to_disease.txt"
echo "== disease names =="
for d in $(awk -F'\t' -v g="$G" '$2==g {print $4}' "$H/genes_to_disease.txt"); do
  grep -m1 -P "^$d\t" "$H/phenotype.hpoa" | cut -f1,2
done
echo "== number of distinct phenotype terms =="
awk -F'\t' -v g="$G" '$2==g {print $3}' "$H/genes_to_phenotype.txt" | sort -u | wc -l
grep -m1 "BS1" /tmp/known_check.acmg.tsv | cut -f22

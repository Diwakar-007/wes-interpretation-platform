#!/usr/bin/env bash
# Shows the gnomAD constraint columns we need (gene, MANE flag, LOEUF) and an example gene.
F="$HOME/ngs-data/annotation/gnomad/gnomad.v4.1.constraint_metrics.tsv"
head -1 "$F" | tr '\t' '\n' | grep -n -E '^(gene|gene_id|transcript|canonical|mane_select)$|^lof\.' | head -20
echo "--- CFTR / BRCA1 rows (gene, transcript, canonical, mane_select) ---"
awk -F'\t' 'NR>1 && ($1=="CFTR" || $1=="BRCA1")' "$F" | cut -f1,3,4,5 | head -6

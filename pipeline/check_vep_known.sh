#!/usr/bin/env bash
# Sanity check: annotate two well-known variants and confirm gnomAD, ClinVar, REVEL/AlphaMissense are filled.
#   MTHFR c.665C>T (rs1801133, common, ~30% allele frequency)
#   CFTR p.Phe508del (rs113993960, the classic cystic fibrosis pathogenic variant)
set -euo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate vep
A="$HOME/ngs-data/annotation"
IN="$(dirname "$0")/assets/known_variants_check.vcf"
OUT=/tmp/known_check.vep.vcf

vep --input_file "$IN" --format vcf --output_file "$OUT" --vcf --force_overwrite --no_stats \
    --offline --cache --dir_cache "$A/vep_cache" --cache_version 116 --assembly GRCh38 \
    --fasta "$HOME/ngs-data/refs/GRCh38/Homo_sapiens_assembly38.fasta" \
    --hgvs --symbol --mane --check_existing --af_gnomade --af_gnomadg --max_af \
    --flag_pick_allele_gene --pick_order mane_select,canonical,appris,tsl,biotype,ccds,rank,length \
    --custom file=$A/clinvar/clinvar.chr.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \
    --dir_plugins "$A/vep_plugins" \
    --plugin AlphaMissense,file=$A/alphamissense/AlphaMissense_hg38.tsv.gz \
    --plugin REVEL,file=$A/revel/revel_grch38_chr.tsv.gz 2>/dev/null

bcftools +split-vep -d -i 'PICK="1"' \
  -f '%ID | %SYMBOL | %HGVSp | %Existing_variation | gnomADe_AF=%gnomADe_AF | gnomADg_AF=%gnomADg_AF | ClinVar=%ClinVar_CLNSIG (%ClinVar_CLNREVSTAT) | REVEL=%REVEL | AM=%am_class\n' \
  "$OUT" 2>/dev/null

#!/usr/bin/env bash
# Downloads GRCh38 (GATK/Broad version) with a ready-made BWA index and BQSR known-sites files.
# Total about 11 GB. It resumes if interrupted: just run it again.
# Usage:  bash scripts/setup/03_download_grch38.sh
#
# Note: for BQSR known SNPs we use 1000G phase1 high-confidence SNPs (1.9 GB) instead of
# full dbSNP (about 11 GB uncompressed). The Broad lists both as valid BQSR known sites.
set -euo pipefail

OUT="$HOME/ngs-data/refs/GRCh38"
B=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
mkdir -p "$OUT"
cd "$OUT"

FILES=(
  # reference genome + index + sequence dictionary
  Homo_sapiens_assembly38.fasta
  Homo_sapiens_assembly38.fasta.fai
  Homo_sapiens_assembly38.dict
  # prebuilt BWA index (the .64 prefix is just Broad's naming; .alt makes bwa ALT-contig aware)
  Homo_sapiens_assembly38.fasta.64.alt
  Homo_sapiens_assembly38.fasta.64.amb
  Homo_sapiens_assembly38.fasta.64.ann
  Homo_sapiens_assembly38.fasta.64.bwt
  Homo_sapiens_assembly38.fasta.64.pac
  Homo_sapiens_assembly38.fasta.64.sa
  # known sites for BQSR
  Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
  Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
  Homo_sapiens_assembly38.known_indels.vcf.gz
  Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
  1000G_phase1.snps.high_confidence.hg38.vcf.gz
  1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
)

for f in "${FILES[@]}"; do
  echo ">> $f"
  wget -c -q --show-progress "$B/$f"
done

# bwa looks for <prefix>.bwt, <prefix>.sa etc. Create links so the prefix is the fasta name.
for ext in alt amb ann bwt pac sa; do
  ln -sf "Homo_sapiens_assembly38.fasta.64.$ext" "Homo_sapiens_assembly38.fasta.$ext"
done

echo
du -sh "$OUT"
echo ">> GRCh38 reference ready in $OUT"

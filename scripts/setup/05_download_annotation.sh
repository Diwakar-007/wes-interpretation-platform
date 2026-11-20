#!/usr/bin/env bash
# Downloads and prepares annotation resources for GRCh38 (all free for academic/non-commercial use).
#
#   cache        Ensembl VEP 116 offline cache (27.6 GB download, ~30 GB extracted; includes gnomAD AFs)
#   clinvar      ClinVar VCF, contigs renamed 1 -> chr1 to match our reference
#   alphamissense AlphaMissense missense pathogenicity scores (0.6 GB)
#   revel        REVEL missense scores, reformatted for GRCh38 + tabix (0.7 GB download, ~6 GB temporary)
#   constraint   gnomAD v4.1 gene constraint (LOEUF) for the ACMG PVS1 rule
#
# Usage:  bash scripts/setup/05_download_annotation.sh            # everything
#         bash scripts/setup/05_download_annotation.sh clinvar    # one resource
# Safe to re-run: finished steps are skipped, interrupted downloads resume.
set -euo pipefail

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate vep        # provides bgzip, tabix, bcftools

RELEASE=116
A="$HOME/ngs-data/annotation"
mkdir -p "$A"
STEPS=${*:-"clinvar constraint alphamissense revel cache"}

done_marker() { [ -f "$1/.done" ]; }

for step in $STEPS; do
case $step in

cache)
  D="$A/vep_cache"; mkdir -p "$D"
  if done_marker "$D"; then echo ">> VEP cache already ready"; continue; fi
  echo ">> VEP cache $RELEASE: downloading 27.6 GB (resumable)..."
  URL="https://ftp.ensembl.org/pub/release-$RELEASE/variation/indexed_vep_cache/homo_sapiens_vep_${RELEASE}_GRCh38.tar.gz"
  if command -v aria2c >/dev/null; then
    # Ensembl limits speed per connection; 16 parallel connections are ~3-10x faster than wget
    aria2c -c -x 16 -s 16 -k 1M --file-allocation=none --summary-interval=60 \
      --max-tries=0 --retry-wait=10 -d "$D" -o "homo_sapiens_vep_${RELEASE}_GRCh38.tar.gz" "$URL"
  else
    wget -c -q --show-progress -P "$D" "$URL"
  fi
  echo ">> Extracting (10-30 minutes)..."
  tar -xzf "$D/homo_sapiens_vep_${RELEASE}_GRCh38.tar.gz" -C "$D"
  rm "$D/homo_sapiens_vep_${RELEASE}_GRCh38.tar.gz"
  touch "$D/.done"
  du -sh "$D"
  ;;

clinvar)
  D="$A/clinvar"; mkdir -p "$D"
  if done_marker "$D"; then echo ">> ClinVar already ready"; continue; fi
  echo ">> ClinVar..."
  wget -q --show-progress -O "$D/clinvar.vcf.gz" https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
  # ClinVar's header has no ##contig lines; indexing lets bcftools learn the contig names
  tabix -f -p vcf "$D/clinvar.vcf.gz"
  # ClinVar uses 1,2,..,MT; our reference uses chr1,chr2,..,chrM
  { for c in $(seq 1 22) X Y; do echo -e "$c\tchr$c"; done; echo -e "MT\tchrM"; } > "$D/chr_rename.txt"
  bcftools annotate --rename-chrs "$D/chr_rename.txt" -Oz -o "$D/clinvar.chr.vcf.gz" "$D/clinvar.vcf.gz"
  tabix -f -p vcf "$D/clinvar.chr.vcf.gz"
  bcftools view -h "$D/clinvar.vcf.gz" | grep -m1 '^##fileDate' > "$D/VERSION.txt"
  rm "$D/clinvar.vcf.gz" "$D/clinvar.vcf.gz.tbi"
  touch "$D/.done"
  echo "   $(cat "$D/VERSION.txt"), $(bcftools index -n "$D/clinvar.chr.vcf.gz") variants"
  ;;

constraint)
  D="$A/gnomad"; mkdir -p "$D"
  if done_marker "$D"; then echo ">> gnomAD constraint already ready"; continue; fi
  echo ">> gnomAD v4.1 constraint..."
  wget -q --show-progress -O "$D/gnomad.v4.1.constraint_metrics.tsv" \
    https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv
  touch "$D/.done"
  ;;

alphamissense)
  D="$A/alphamissense"; mkdir -p "$D"
  if done_marker "$D"; then echo ">> AlphaMissense already ready"; continue; fi
  echo ">> AlphaMissense..."
  wget -c -q --show-progress -O "$D/AlphaMissense_hg38.tsv.gz" \
    https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
  # The file has 3 copyright lines + 1 header line; tabix needs bgzip compression.
  if ! tabix -f -s 1 -b 2 -e 2 -S 4 "$D/AlphaMissense_hg38.tsv.gz" 2>/dev/null; then
    echo "   re-compressing with bgzip..."
    zcat "$D/AlphaMissense_hg38.tsv.gz" | bgzip -@ 4 > "$D/AlphaMissense_hg38.bgz.tsv.gz"
    mv "$D/AlphaMissense_hg38.bgz.tsv.gz" "$D/AlphaMissense_hg38.tsv.gz"
    tabix -f -s 1 -b 2 -e 2 -S 4 "$D/AlphaMissense_hg38.tsv.gz"
  fi
  touch "$D/.done"
  ;;

revel)
  D="$A/revel"; mkdir -p "$D/tmp"
  if done_marker "$D"; then echo ">> REVEL already ready"; continue; fi
  echo ">> REVEL: downloading..."
  wget -c -q --show-progress -O "$D/revel-v1.3_all_chromosomes.zip" \
    https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip
  echo ">> REVEL: converting to GRCh38, tab-separated, chr-prefixed, sorted (15-30 minutes)..."
  # columns: chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid
  # '|| true': head closes the pipe early, which makes unzip exit non-zero (harmless here)
  HEADER=$(unzip -p "$D/revel-v1.3_all_chromosomes.zip" | head -1 | tr ',' '\t' | sed 's/^/#/' || true)
  {
    echo "$HEADER"
    unzip -p "$D/revel-v1.3_all_chromosomes.zip" | tail -n +2 | tr ',' '\t' \
      | awk 'BEGIN{FS=OFS="\t"} $3!="." {$1="chr"$1; print}' \
      | sort -t $'\t' -k1,1 -k3,3n -S 2G --parallel=4 -T "$D/tmp"
  } | bgzip -@ 4 > "$D/revel_grch38_chr.tsv.gz"
  tabix -f -s 1 -b 3 -e 3 "$D/revel_grch38_chr.tsv.gz"
  rm -rf "$D/tmp" "$D/revel-v1.3_all_chromosomes.zip"
  touch "$D/.done"
  ;;

hpo)
  # Human Phenotype Ontology: ontology + gene/disease annotations (incl. inheritance modes)
  D="$A/hpo"; mkdir -p "$D"
  if done_marker "$D"; then echo ">> HPO already ready"; continue; fi
  echo ">> HPO..."
  H=https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download
  for f in hp.obo genes_to_phenotype.txt genes_to_disease.txt phenotype.hpoa; do
    wget -q -O "$D/$f" "$H/$f"
  done
  grep -m1 '^data-version' "$D/hp.obo" > "$D/VERSION.txt"
  touch "$D/.done"
  echo "   $(cat "$D/VERSION.txt")"
  ;;

*) echo "Unknown step: $step (use: cache clinvar constraint alphamissense revel hpo)"; exit 1 ;;
esac
done

echo
du -sh "$A"/* 2>/dev/null
echo ">> Annotation resources in $A"

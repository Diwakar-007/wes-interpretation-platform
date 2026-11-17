#!/usr/bin/env bash
# Downloads a tiny human test dataset (a small piece of chr22, a few MB) from nf-core.
# We use it to build and debug the pipeline in minutes instead of hours.
# Usage:  conda activate ngs && bash scripts/setup/02_download_test_data.sh
set -euo pipefail

OUT="$HOME/ngs-data/test"
T=https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens
mkdir -p "$OUT"/{ref,fastq}

echo ">> Downloading test reference + known sites..."
for f in genome.fasta genome.fasta.fai genome.dict genome.bed; do
  wget -q -c -O "$OUT/ref/$f" "$T/genome/$f"
done
for f in dbsnp_146.hg38.vcf.gz dbsnp_146.hg38.vcf.gz.tbi mills_and_1000G.indels.vcf.gz mills_and_1000G.indels.vcf.gz.tbi; do
  wget -q -c -O "$OUT/ref/$f" "$T/genome/vcf/$f"
done

echo ">> Downloading test FASTQ pair..."
wget -q -c -O "$OUT/fastq/test_R1.fastq.gz" "$T/illumina/fastq/test_1.fastq.gz"
wget -q -c -O "$OUT/fastq/test_R2.fastq.gz" "$T/illumina/fastq/test_2.fastq.gz"

echo ">> Building BWA index for the test genome (seconds)..."
bwa index "$OUT/ref/genome.fasta" 2> "$OUT/ref/bwa_index.log"

echo
ls -lh "$OUT"/ref "$OUT"/fastq
echo ">> Test data ready in $OUT"

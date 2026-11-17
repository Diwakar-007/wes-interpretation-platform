#!/usr/bin/env bash
# Checks that the reference and test-data URLs exist and prints their sizes.
B=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
for f in Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.dict \
         Homo_sapiens_assembly38.fasta.64.alt Homo_sapiens_assembly38.fasta.64.amb Homo_sapiens_assembly38.fasta.64.ann \
         Homo_sapiens_assembly38.fasta.64.bwt Homo_sapiens_assembly38.fasta.64.pac Homo_sapiens_assembly38.fasta.64.sa \
         Mills_and_1000G_gold_standard.indels.hg38.vcf.gz Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
         Homo_sapiens_assembly38.known_indels.vcf.gz Homo_sapiens_assembly38.known_indels.vcf.gz.tbi \
         1000G_phase1.snps.high_confidence.hg38.vcf.gz 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi; do
  s=$(curl -sIL "$B/$f" | grep -i '^content-length' | tail -1 | tr -dc 0-9)
  echo "GRCh38  ${s:-0} bytes  $f"
done

T=https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens
for f in genome/genome.fasta genome/genome.fasta.fai genome/genome.dict genome/genome.bed \
         genome/vcf/dbsnp_146.hg38.vcf.gz genome/vcf/dbsnp_146.hg38.vcf.gz.tbi \
         genome/vcf/mills_and_1000G.indels.vcf.gz genome/vcf/mills_and_1000G.indels.vcf.gz.tbi \
         illumina/fastq/test_1.fastq.gz illumina/fastq/test_2.fastq.gz; do
  c=$(curl -sIL -o /dev/null -w '%{http_code}' "$T/$f")
  echo "TEST    HTTP $c  $f"
done

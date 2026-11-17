#!/usr/bin/env bash
# Shows mapping stats and the top contigs for a BAM.  Usage: bash debug_bam.sh <bam>
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ngs
BAM=$1
samtools flagstat "$BAM" | head -8
echo "== top contigs by mapped reads =="
samtools idxstats "$BAM" 2>/dev/null | sort -k3,3nr | head -8 || \
  samtools view -F 4 "$BAM" | cut -f3 | sort | uniq -c | sort -nr | head -8
echo "== read group in header =="
samtools view -H "$BAM" | grep '^@RG'

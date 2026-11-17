#!/usr/bin/env bash
# Prints tool versions, data status and disk/RAM. Handy to paste when asking for help.
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ngs
export PATH="$HOME/.local/bin:$PATH"
echo "== tools =="
echo "nextflow : $(nextflow -v 2>/dev/null)"
echo "bwa      : $(bwa 2>&1 | grep -i '^version')"
echo "samtools : $(samtools --version | head -1)"
echo "bcftools : $(bcftools --version | head -1)"
echo "gatk     : $(gatk --version 2>/dev/null | grep -m1 -i gatk)"
echo "fastp    : $(fastp --version 2>&1)"
echo "mosdepth : $(mosdepth --version)"
echo "docker   : $(docker --version)"
echo "== test data =="
ls -lh ~/ngs-data/test/ref ~/ngs-data/test/fastq 2>&1
echo "== GRCh38 =="
ls -lh ~/ngs-data/refs/GRCh38 | head -30
echo "== resources =="
df -h ~ | tail -1
free -h | head -2

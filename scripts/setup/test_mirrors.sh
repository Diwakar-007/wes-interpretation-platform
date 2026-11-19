#!/usr/bin/env bash
# Tests download speed of the VEP cache from different servers (10 seconds each, nothing saved).
F=release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz
for u in "https://ftp.ensembl.org/pub/$F" "https://ftp.ebi.ac.uk/ensemblorg/pub/$F" "http://ftp.ensembl.org/pub/$F"; do
  sp=$(curl -s -o /dev/null -m 10 -w '%{speed_download}' "$u")
  printf "%6.2f MB/s  %s\n" "$(echo "$sp/1048576" | bc -l)" "$u"
done

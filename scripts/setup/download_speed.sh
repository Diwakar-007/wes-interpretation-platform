#!/usr/bin/env bash
# Measures how fast the annotation downloads are growing over 15 seconds.
A="$HOME/ngs-data/annotation"
snap() { find "$A" -type f \( -name '*.tar.gz' -o -name '*.zip' -o -name '*.tsv*' -o -name '*.vcf.gz' \) -printf '%s %p\n' 2>/dev/null | sort -k2; }
before=$(snap); sleep 15; after=$(snap)
echo "file | size now | speed"
join -1 2 -2 2 <(echo "$before") <(echo "$after") | while read -r p s1 s2; do
  printf "%s | %d MB | %.1f MB/s\n" "${p#$A/}" $((s2/1048576)) "$(echo "($s2-$s1)/15/1048576" | bc -l)"
done
TOTAL=27644657162
CUR=$(stat -c %s "$A/vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz" 2>/dev/null || echo 0)
echo "VEP cache: $((CUR*100/TOTAL))% of 27.6 GB"
echo "steps finished: $(grep -a '>>' "$A/download.log" 2>/dev/null | tr '\n' ' ')"

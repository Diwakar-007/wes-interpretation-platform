#!/usr/bin/env bash
# Finds the current Ensembl release, the matching bioconda VEP version and download sizes.
size() { curl -sIL "$1" | grep -i '^content-length' | tail -1 | tr -dc 0-9; }

REL=$(curl -sL https://ftp.ensembl.org/pub/current_variation/indexed_vep_cache/ \
      | grep -o 'homo_sapiens_vep_[0-9]*_GRCh38.tar.gz' | head -1 | sed -E 's/homo_sapiens_vep_([0-9]+)_GRCh38.tar.gz/\1/')
echo "Ensembl release: $REL"
echo "VEP cache size : $(size https://ftp.ensembl.org/pub/release-$REL/variation/indexed_vep_cache/homo_sapiens_vep_${REL}_GRCh38.tar.gz) bytes"

echo "== bioconda ensembl-vep versions (latest 5) =="
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda search -c conda-forge -c bioconda ensembl-vep 2>/dev/null | awk 'NR>2{print $2}' | sort -uV | tail -5

echo "== plugins (release/$REL branch) =="
for p in AlphaMissense REVEL; do
  echo "$(curl -s -o /dev/null -w '%{http_code}' https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/$REL/$p.pm)  $p.pm"
done

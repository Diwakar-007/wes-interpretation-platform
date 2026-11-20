#!/usr/bin/env bash
# Shows VEP cache download / extraction progress.
D="$HOME/ngs-data/annotation/vep_cache"
F="$D/homo_sapiens_vep_116_GRCh38.tar.gz"
TOTAL=27644657162
if [ -f "$D/.done" ]; then echo "CACHE DONE: $(du -sh "$D" | cut -f1)"; exit 0; fi
echo "== processes =="
pgrep -a aria2c | cut -c1-60 || echo "aria2c: not running"
pgrep -a wget | cut -c1-60 || true
pgrep -a -x tar || echo "tar (extracting): not running"
echo "== files =="
ls -la "$D"
if [ -f "$F" ]; then
  # with --file-allocation=none the file size grows; the .aria2 file means download not finished
  S1=$(du -b "$F" | cut -f1); sleep 20; S2=$(du -b "$F" | cut -f1)
  echo "downloaded: $((S2/1048576)) MB of $((TOTAL/1048576)) MB ($((S2*100/TOTAL))%)"
  SPEED=$(( (S2-S1)/20 ))
  echo "speed: $((SPEED/1024)) KB/s"
  [ "$SPEED" -gt 0 ] && echo "ETA: $(( (TOTAL-S2)/SPEED/60 )) minutes"
  [ -f "$F.aria2" ] && echo "(.aria2 control file present: download in progress or paused)"
fi
[ -d "$D/homo_sapiens" ] && echo "extracted so far: $(du -sh "$D/homo_sapiens" | cut -f1)"

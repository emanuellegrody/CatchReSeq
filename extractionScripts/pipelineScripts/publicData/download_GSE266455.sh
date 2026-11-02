#!/bin/bash
# Download published GSE266455 (HIV-Seq, GRCh38 + HIV clade B consensus) 10x matrices.
# One subdirectory per sample holding standard 10x filenames, so Read10X() works directly.

set -euo pipefail

BASEURL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266455/suppl"
OUTDIR="/projects/b1042/GoyalLab/egrody/publicData/HIVseq/matrices"
LOGDIR="/projects/b1042/GoyalLab/egrody/publicData/HIVseq/scripts/logs"
mkdir -p "$OUTDIR" "$LOGDIR"

SAMPLES=(
  PID1052W0w1NOTHING
  PID1052W0w3HIV
  PID1052W0w5HIV
  PID1052W24w2NOTHING
  PID1052W24w4HIV
  PID1052W24w6HIV
  PID145W0w1HIV
  PID145W0w2HIV
  PID8026W0w5HIV
  PID8026W0w6HIV
  PID8026W0w7HIV
  PID8026W0w8HIV
  PID8026W24w1HIV
  PID8026W24w2HIV
  PID8026W24w3HIV
  PID8026W24w4HIV
  PID8027W0w1NOTHING
  PID8027W0w2NOTHING
  PID8027W0w5HIV
  PID8027W0w6HIV
  PID8027W45w3NOTHING
  PID8027W45w4NOTHING
  PID8027W45w7HIV
  PID8027W45w8HIV
)

for sample in "${SAMPLES[@]}"; do
  sampledir="${OUTDIR}/${sample}"
  mkdir -p "$sampledir"
  for part in barcodes features matrix; do
    case "$part" in
      matrix) target="${sampledir}/matrix.mtx.gz" ;;
      *)      target="${sampledir}/${part}.tsv.gz" ;;
    esac
    if [[ -s "$target" ]]; then
      echo "[skip] $sample $part"
      continue
    fi
    case "$part" in
      matrix) remote="GSE266455_SCOPE_Samples_${sample}_matrix.mtx.gz" ;;
      *)      remote="GSE266455_SCOPE_Samples_${sample}_${part}.tsv.gz" ;;
    esac
    echo "[get ] $sample $part"
    wget -q --tries=5 --waitretry=10 -O "${target}.part" "${BASEURL}/${remote}"
    mv "${target}.part" "$target"
  done
done

echo "=== verifying gzip integrity ==="
find "$OUTDIR" -name "*.gz" -print0 | xargs -0 -n1 -P4 gzip -t
echo "=== done: $(find "$OUTDIR" -name "*.gz" | wc -l) files ==="
du -sh "$OUTDIR"

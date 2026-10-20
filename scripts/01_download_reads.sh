#!/usr/bin/env bash
# =====================================================================
# 01  Download raw reads (paired-end FASTQ) for the 8 airway samples from SRA
#
# Usage:  bash scripts/01_download_reads.sh [READS_PER_SAMPLE]
#   READS_PER_SAMPLE  read pairs to keep per sample (default 2000000).
#                     Use 0 to download everything (~20-30 million pairs per
#                     sample, roughly 30 GB in total).
#   2 million pairs per sample is plenty to reproduce the main result on a laptop.
#
# Output: data/fastq/<run>_1.fastq.gz and <run>_2.fastq.gz
# =====================================================================
set -euo pipefail

READS=${1:-2000000}
OUT=data/fastq
mkdir -p "$OUT"

# Run IDs from the sample table (first column, skipping the header)
RUNS=$(tail -n +2 data/samples.tsv | cut -f1)

for run in $RUNS; do
  if [[ -s "$OUT/${run}_1.fastq.gz" && -s "$OUT/${run}_2.fastq.gz" ]]; then
    echo "[skip] $run already downloaded"
    continue
  fi
  echo "[download] $run"
  if [[ "$READS" -gt 0 ]]; then
    # -X: only the first N spots (read pairs); --split-files: separate R1 / R2
    fastq-dump --split-files --gzip -X "$READS" --outdir "$OUT" "$run"
  else
    prefetch "$run" --output-directory "$OUT/sra"
    fasterq-dump --split-files --threads 4 --outdir "$OUT" "$OUT/sra/$run"
    pigz -p 4 "$OUT/${run}"_*.fastq
    rm -rf "$OUT/sra/$run"
  fi
done

echo "Done. Files in $OUT:"
ls -lh "$OUT"/*.fastq.gz

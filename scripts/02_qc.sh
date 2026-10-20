#!/usr/bin/env bash
# =====================================================================
# 02  Read quality control: FastQC on every file, summarised by MultiQC
#
# Usage:  bash scripts/02_qc.sh [THREADS]
# Output: results/qc/fastqc/*.html  and  results/qc/multiqc_report.html
#
# What to look at in the MultiQC report:
#   - per-base quality (should stay high along the read)
#   - adapter content (low for this library type)
#   - duplication levels (RNA-seq is naturally high for very abundant genes)
#   - GC content (one smooth peak; a second peak can mean contamination)
# =====================================================================
set -euo pipefail

THREADS=${1:-4}
mkdir -p results/qc/fastqc

fastqc --threads "$THREADS" --outdir results/qc/fastqc data/fastq/*.fastq.gz
multiqc --force --outdir results/qc results/qc/fastqc

echo "Open results/qc/multiqc_report.html in a web browser."

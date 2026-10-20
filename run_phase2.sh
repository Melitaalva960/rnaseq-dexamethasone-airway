#!/usr/bin/env bash
# Phase 2 in one go: SRA reads -> QC -> Salmon -> tximport -> DESeq2 -> figures -> comparison.
# Run from the project folder in Linux / WSL, after:
#   conda env create -f envs/phase2.yml && conda activate rnaseq-airway
# and after Phase 1 (run_phase1.R) so there is something to compare with.
#
# Usage:  bash run_phase2.sh [READS_PER_SAMPLE] [THREADS]
set -euo pipefail
READS=${1:-2000000}
THREADS=${2:-4}

bash scripts/01_download_reads.sh "$READS"
bash scripts/02_qc.sh "$THREADS"
bash scripts/03_salmon_quant.sh "$THREADS"

Rscript R/01_get_counts.R input=salmon out_dir=results/phase2/counts
Rscript R/02_deseq2.R input=salmon counts_dir=results/phase2/counts out_dir=results/phase2
Rscript R/03_figures.R dds=results/phase2/dds.rds table=results/phase2/tables/deseq2_all_genes.csv \
  out_dir=results/phase2/figures
Rscript R/04_enrichment.R table=results/phase2/tables/deseq2_all_genes.csv \
  table_dir=results/phase2/tables fig_dir=results/phase2/figures
Rscript R/05_compare_phases.R

echo "Done. Phase 2 results are in results/phase2; comparison in results/figures/phase_comparison.png"

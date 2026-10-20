#!/usr/bin/env bash
# =====================================================================
# 03  Transcript quantification with Salmon
#
# 1. Downloads the GENCODE human transcriptome (once)
# 2. Builds a Salmon index (once)
# 3. Quantifies each sample (selective alignment, automatic library type)
# 4. Writes a transcript -> gene table for tximport
#
# Usage:  bash scripts/03_salmon_quant.sh [THREADS]
#
# Settings (environment variables):
#   GENCODE_RELEASE=46        GENCODE human release to use
#   PROTEIN_CODING_ONLY=1     index only protein-coding transcripts. Uses much
#                             less memory (good for laptops with 8 GB RAM);
#                             set to 0 for the full transcriptome (~16 GB RAM).
#
# Output: results/salmon/<run>/quant.sf, data/reference/tx2gene.tsv
# =====================================================================
set -euo pipefail

THREADS=${1:-4}
GENCODE_RELEASE=${GENCODE_RELEASE:-46}
PROTEIN_CODING_ONLY=${PROTEIN_CODING_ONLY:-1}

REF=data/reference
mkdir -p "$REF" results/salmon

FASTA_GZ="$REF/gencode.v${GENCODE_RELEASE}.transcripts.fa.gz"
URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}/gencode.v${GENCODE_RELEASE}.transcripts.fa.gz"

# ---- 1. Transcriptome --------------------------------------------------------
if [[ ! -s "$FASTA_GZ" ]]; then
  echo "[download] GENCODE v${GENCODE_RELEASE} transcripts"
  curl -L --fail -o "$FASTA_GZ" "$URL"
fi

TX_FA="$REF/transcripts.fa"
if [[ ! -s "$TX_FA" ]]; then
  if [[ "$PROTEIN_CODING_ONLY" == "1" ]]; then
    echo "[filter] keeping protein-coding transcripts only"
    # GENCODE headers end with |<transcript_biotype>|
    zcat "$FASTA_GZ" | awk '/^>/ {keep = ($0 ~ /\|protein_coding\|$/)} keep' > "$TX_FA"
  else
    zcat "$FASTA_GZ" > "$TX_FA"
  fi
fi

# ---- 2. transcript -> gene table (from the GENCODE FASTA headers) ------------
# Header: >ENST...|ENSG...|OTTHUMG...|OTTHUMT...|transcript_name|gene_name|length|biotype|
if [[ ! -s "$REF/tx2gene.tsv" ]]; then
  {
    printf "transcript_id\tgene_id\tgene_name\n"
    # version numbers removed (ENST...1.2 -> ENST...1) so IDs match across tools
    grep '^>' "$TX_FA" | sed 's/^>//' | awk -F'|' -v OFS='\t' '{sub(/\.[0-9]+$/, "", $1); sub(/\.[0-9]+$/, "", $2); print $1, $2, $6}'
  } > "$REF/tx2gene.tsv"
  echo "[tx2gene] $(($(wc -l < "$REF/tx2gene.tsv") - 1)) transcripts"
fi

# ---- 3. Salmon index ---------------------------------------------------------------
INDEX="$REF/salmon_index"
if [[ ! -d "$INDEX" ]]; then
  echo "[index] building Salmon index (can take 10-40 minutes)"
  salmon index --transcripts "$TX_FA" --index "$INDEX" --gencode --threads "$THREADS"
fi

# ---- 4. Quantify each sample ----------------------------------------------------------
for run in $(tail -n +2 data/samples.tsv | cut -f1); do
  if [[ -s "results/salmon/$run/quant.sf" ]]; then
    echo "[skip] $run already quantified"; continue
  fi
  echo "[quant] $run"
  salmon quant --index "$INDEX" --libType A \
    --mates1 "data/fastq/${run}_1.fastq.gz" --mates2 "data/fastq/${run}_2.fastq.gz" \
    --validateMappings --gcBias --threads "$THREADS" \
    --output "results/salmon/$run"
done

# Mapping rate per sample (also shown in MultiQC if you re-run it on results/)
echo; echo "Mapping rates:"
for run in $(tail -n +2 data/samples.tsv | cut -f1); do
  rate=$(grep -o '"percent_mapped": [0-9.]*' "results/salmon/$run/aux_info/meta_info.json" | awk '{print $2}')
  printf "  %s  %.1f%%\n" "$run" "$rate"
done

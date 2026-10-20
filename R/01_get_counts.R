## =====================================================================
## 01  Gene counts and sample table
##
##  input=airway  (Phase 1) counts from the Bioconductor 'airway' package
##                (Himes et al. 2014, GEO GSE52778; reads aligned and
##                counted per Ensembl gene by the package authors)
##  input=salmon  (Phase 2) your own Salmon quantification, summarised to
##                genes with tximport
##
## Output: results/counts/gene_counts.tsv   genes x samples (integer counts)
##         results/counts/sample_info.tsv   one row per sample
##         results/counts/txi.rds           (salmon only) tximport object
## =====================================================================
source("R/utils.R")

s <- get_settings(list(
  input      = "airway",
  salmon_dir = "results/salmon",           # one folder per run, each with quant.sf
  tx2gene    = "data/reference/tx2gene.tsv",
  samples    = "data/samples.tsv",
  out_dir    = "results/counts"
))
ensure_dir(s$out_dir)

## Sample table (run ID, cell line, treatment) shared by both phases
samples <- read.delim(s$samples)
samples$treatment <- factor(samples$treatment, levels = c("untreated", "dexamethasone"))

if (s$input == "airway") {
  need(c("airway", "SummarizedExperiment"))
  data("airway", package = "airway")
  counts <- SummarizedExperiment::assay(airway, "counts")
  stopifnot("Runs in data/samples.tsv must match the airway package" =
              setequal(colnames(counts), samples$run))
  counts <- counts[, samples$run]

} else if (s$input == "salmon") {
  need("tximport")
  files <- file.path(s$salmon_dir, samples$run, "quant.sf")
  names(files) <- samples$run
  if (!all(file.exists(files))) {
    stop("Missing Salmon output for: ", paste(samples$run[!file.exists(files)], collapse = ", "),
         "\nRun scripts/03_salmon_quant.sh first.")
  }
  tx2gene <- read.delim(s$tx2gene, header = TRUE)[, c("transcript_id", "gene_id")]
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
  saveRDS(txi, file.path(s$out_dir, "txi.rds"))
  counts <- round(txi$counts)

} else {
  stop("input must be 'airway' or 'salmon'")
}

## Ensembl gene IDs without version suffix (ENSG00000123456.7 -> ENSG00000123456)
rownames(counts) <- sub("\\.\\d+$", "", rownames(counts))
storage.mode(counts) <- "integer"

write.table(data.frame(gene_id = rownames(counts), counts, check.names = FALSE),
            file.path(s$out_dir, "gene_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(samples, file.path(s$out_dir, "sample_info.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("%s counts: %d genes x %d samples; %.1f million reads per sample on average\n",
            s$input, nrow(counts), ncol(counts), mean(colSums(counts)) / 1e6))
cat("Saved to", s$out_dir, "\n")

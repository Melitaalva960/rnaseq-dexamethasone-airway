## =====================================================================
## 02  Differential expression with DESeq2
##
## Design: ~ cell_line + treatment
##   Each cell line was measured treated and untreated, so cell line is
##   included as a blocking factor (a paired analysis). This removes
##   donor-to-donor differences and gives more power to detect the drug effect.
##
## Output: results/tables/deseq2_all_genes.csv      every tested gene
##         results/tables/deseq2_significant.csv    padj < 0.05 and |log2FC| >= 1
##         results/tables/normalized_counts.csv
##         results/tables/summary.txt
##         results/dds.rds                          for the figure scripts
## =====================================================================
source("R/utils.R")
need(c("DESeq2", "AnnotationDbi", "org.Hs.eg.db"))

s <- get_settings(list(
  counts_dir = "results/counts",
  input      = "airway",          # "salmon" uses txi.rds so transcript lengths are accounted for
  out_dir    = "results",
  min_count  = 10,                # keep genes with >= min_count reads ...
  min_samples = 4                 # ... in at least this many samples (size of the smallest group)
))
tables <- file.path(s$out_dir, "tables"); ensure_dir(tables)

samples <- read.delim(file.path(s$counts_dir, "sample_info.tsv"))
samples$treatment <- factor(samples$treatment, levels = c("untreated", "dexamethasone"))
samples$cell_line <- factor(samples$cell_line)
rownames(samples) <- samples$run

if (s$input == "salmon") {
  txi <- readRDS(file.path(s$counts_dir, "txi.rds"))
  txi <- lapply(txi, function(x) if (is.matrix(x)) x[, samples$run] else x)
  dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ cell_line + treatment)
  rownames(dds) <- sub("\\.\\d+$", "", rownames(dds))
} else {
  counts <- as.matrix(read.delim(file.path(s$counts_dir, "gene_counts.tsv"), row.names = 1, check.names = FALSE))
  dds <- DESeqDataSetFromMatrix(counts[, samples$run], colData = samples, design = ~ cell_line + treatment)
}

## Pre-filter genes that are barely expressed: speeds things up and
## reduces the multiple-testing burden
n_before <- nrow(dds)
dds <- dds[rowSums(counts(dds) >= s$min_count) >= s$min_samples, ]
cat(sprintf("Kept %d of %d genes after filtering\n", nrow(dds), n_before))

dds <- DESeq(dds)

## Raw results (for p-values) and shrunken log2 fold changes (for ranking
## and plotting: shrinkage stops noisy low-count genes having huge fold changes)
res <- results(dds, name = "treatment_dexamethasone_vs_untreated", alpha = PADJ_MAX)
shrink_type <- if (requireNamespace("apeglm", quietly = TRUE)) "apeglm" else "normal"
if (shrink_type != "apeglm") message("apeglm not installed, using type = 'normal' shrinkage")
res_shr <- lfcShrink(dds, coef = "treatment_dexamethasone_vs_untreated", type = shrink_type, quiet = TRUE)

## Gene symbols and descriptions
ids <- rownames(res)
symbol <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, ids, "SYMBOL", "ENSEMBL", multiVals = "first"))
gname  <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, ids, "GENENAME", "ENSEMBL", multiVals = "first"))

out <- data.frame(
  gene_id = ids, symbol = unname(symbol), description = unname(gname),
  baseMean = res$baseMean,
  log2FC = res_shr$log2FoldChange,          # shrunken
  log2FC_unshrunk = res$log2FoldChange,
  lfcSE = res$lfcSE, stat = res$stat,
  pvalue = res$pvalue, padj = res$padj
)
out <- out[order(out$padj, -abs(out$log2FC), na.last = TRUE), ]
out$direction <- ifelse(is.na(out$padj) | out$padj >= PADJ_MAX | abs(out$log2FC) < LFC_MIN, "not significant",
                        ifelse(out$log2FC > 0, "up", "down"))

sig <- out[out$direction != "not significant", ]
write.csv(out, file.path(tables, "deseq2_all_genes.csv"), row.names = FALSE)
write.csv(sig, file.path(tables, "deseq2_significant.csv"), row.names = FALSE)

norm <- counts(dds, normalized = TRUE)
write.csv(data.frame(gene_id = rownames(norm), symbol = unname(symbol[rownames(norm)]), round(norm, 2),
                     check.names = FALSE),
          file.path(tables, "normalized_counts.csv"), row.names = FALSE)
saveRDS(dds, file.path(s$out_dir, "dds.rds"))

## ---- Summary --------------------------------------------------------------
markers <- out[match(MARKER_GENES, out$symbol), c("symbol", "log2FC", "padj")]
lines <- c(
  sprintf("Input: %s | design: ~ cell_line + treatment | LFC shrinkage: %s", s$input, shrink_type),
  sprintf("Genes tested: %d", nrow(out)),
  sprintf("padj < %.2f (any fold change): %d", PADJ_MAX, sum(out$padj < PADJ_MAX, na.rm = TRUE)),
  sprintf("padj < %.2f and |log2FC| >= %g: %d  (up %d, down %d)", PADJ_MAX, LFC_MIN,
          nrow(sig), sum(sig$direction == "up"), sum(sig$direction == "down")),
  "", "Top 10 genes by adjusted p-value:",
  capture.output(print(head(out[, c("symbol", "log2FC", "padj")], 10), row.names = FALSE, digits = 3)),
  "", "Known glucocorticoid-responsive genes (positive control):",
  capture.output(print(markers, row.names = FALSE, digits = 3))
)
writeLines(lines, file.path(tables, "summary.txt"))
cat(lines, sep = "\n")

## =====================================================================
## 03  Figures
##
##  pca.png                 samples on the first two principal components
##  sample_distances.png    how similar the samples are to each other
##  ma_plot.png             fold change vs mean expression
##  volcano.png             fold change vs significance, top genes labelled
##  top50_heatmap.png       top 50 genes, z-scored expression
##  marker_genes.png        known glucocorticoid-responsive genes, paired by cell line
## =====================================================================
source("R/utils.R")
need(c("DESeq2", "ggplot2", "ggrepel", "pheatmap"))

s <- get_settings(list(
  dds     = "results/dds.rds",
  table   = "results/tables/deseq2_all_genes.csv",
  out_dir = "results/figures"
))
ensure_dir(s$out_dir)
save_plot <- function(p, name, w = 7, h = 5.5) {
  ggsave(file.path(s$out_dir, name), p, width = w, height = h, dpi = 300, bg = "white")
}

dds <- readRDS(s$dds)
res <- read.csv(s$table)
## Variance-stabilised values for visualisation (blind to the design).
## vst() is the fast version and needs >= 1000 genes; fall back for small test sets.
vsd <- if (nrow(dds) >= 1000) vst(dds, blind = TRUE) else varianceStabilizingTransformation(dds, blind = TRUE)

## ---- PCA ---------------------------------------------------------------------
pca <- plotPCA(vsd, intgroup = c("treatment", "cell_line"), returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, colour = treatment, shape = cell_line)) +
  geom_point(size = 4) +
  scale_colour_manual(values = TREAT_COLOURS) +
  labs(title = "PCA of the 500 most variable genes",
       subtitle = "Colour = treatment, shape = cell line",
       x = sprintf("PC1 (%d%% of variance)", pv[1]), y = sprintf("PC2 (%d%% of variance)", pv[2]),
       colour = "Treatment", shape = "Cell line") +
  theme_report()
save_plot(p, "pca.png")

## ---- Sample-to-sample distances ---------------------------------------------------
d <- as.matrix(dist(t(SummarizedExperiment::assay(vsd))))
lab <- paste(dds$cell_line, ifelse(dds$treatment == "dexamethasone", "dex", "ctrl"))
dimnames(d) <- list(lab, lab)
ann <- data.frame(Treatment = dds$treatment, row.names = lab)
pheatmap::pheatmap(d, color = colorRampPalette(c("#08306b", "#6baed6", "white"))(100),
                   annotation_row = ann, annotation_col = ann,
                   annotation_colors = list(Treatment = TREAT_COLOURS),
                   main = "Sample-to-sample distances (VST)",
                   filename = file.path(s$out_dir, "sample_distances.png"), width = 7.5, height = 6)

## ---- MA plot -------------------------------------------------------------------
res$sig <- factor(res$direction, levels = c("up", "down", "not significant"))
dir_cols <- c(up = "#C2410C", down = "#1D4ED8", `not significant` = "grey75")
p <- ggplot(res[order(res$sig, decreasing = TRUE), ], aes(log10(baseMean), log2FC, colour = sig)) +
  geom_point(size = 0.6, alpha = 0.7) +
  geom_hline(yintercept = 0, colour = "grey30") +
  scale_colour_manual(values = dir_cols, name = NULL) +
  labs(title = "MA plot", subtitle = sprintf("Coloured: padj < %.2f and |log2FC| >= %g", PADJ_MAX, LFC_MIN),
       x = "log10 mean normalised count", y = "log2 fold change (dex vs untreated, shrunken)") +
  theme_report()
save_plot(p, "ma_plot.png")

## ---- Volcano -------------------------------------------------------------------------
vol <- res[!is.na(res$padj), ]
vol$neglog10 <- -log10(pmax(vol$padj, 1e-300))
lab_genes <- unique(c(head(vol$symbol[vol$direction == "up"], 10),
                      head(vol$symbol[vol$direction == "down"], 10),
                      intersect(MARKER_GENES, vol$symbol)))
vol$label <- ifelse(vol$symbol %in% lab_genes, vol$symbol, NA)
n_up <- sum(vol$direction == "up"); n_down <- sum(vol$direction == "down")
p <- ggplot(vol[order(vol$sig, decreasing = TRUE), ], aes(log2FC, neglog10, colour = sig)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_vline(xintercept = c(-LFC_MIN, LFC_MIN), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(PADJ_MAX), linetype = "dashed", colour = "grey50") +
  ggrepel::geom_text_repel(aes(label = label), size = 3, colour = "black", max.overlaps = 30,
                           min.segment.length = 0, na.rm = TRUE) +
  scale_colour_manual(values = dir_cols, name = NULL) +
  labs(title = "Dexamethasone vs untreated",
       subtitle = sprintf("%d genes up, %d down (padj < %.2f, |log2FC| >= %g)", n_up, n_down, PADJ_MAX, LFC_MIN),
       x = "log2 fold change (shrunken)", y = "-log10 adjusted p-value") +
  theme_report()
save_plot(p, "volcano.png", 8, 6.5)

## ---- Top 50 heatmap -----------------------------------------------------------------
top <- head(res[res$direction != "not significant", ], 50)
m <- SummarizedExperiment::assay(vsd)[top$gene_id, ]
m <- t(scale(t(m)))                                   # z-score per gene
rownames(m) <- ifelse(is.na(top$symbol), top$gene_id, top$symbol)
colnames(m) <- lab
ord <- order(dds$treatment, dds$cell_line)
ann_col <- data.frame(Treatment = dds$treatment, `Cell line` = dds$cell_line,
                      row.names = lab, check.names = FALSE)
pheatmap::pheatmap(m[, ord], cluster_cols = FALSE,
                   color = colorRampPalette(c("#1D4ED8", "white", "#C2410C"))(101),
                   breaks = seq(-2.5, 2.5, length.out = 102),
                   annotation_col = ann_col, annotation_colors = list(Treatment = TREAT_COLOURS),
                   gaps_col = sum(dds$treatment == "untreated"), fontsize_row = 7,
                   main = "Top 50 differentially expressed genes (z-score of VST expression)",
                   filename = file.path(s$out_dir, "top50_heatmap.png"), width = 7.5, height = 10)

## ---- Marker genes: paired dot plot --------------------------------------------------
gene_ids <- res$gene_id[match(MARKER_GENES, res$symbol)]
keep <- !is.na(gene_ids)
nc <- counts(dds, normalized = TRUE)[gene_ids[keep], , drop = FALSE]
long <- data.frame(
  gene = rep(MARKER_GENES[keep], times = ncol(nc)),
  count = as.vector(nc),
  treatment = rep(dds$treatment, each = nrow(nc)),
  cell_line = rep(dds$cell_line, each = nrow(nc))
)
long$gene <- factor(long$gene, levels = MARKER_GENES[keep])
p <- ggplot(long, aes(treatment, count + 1)) +
  geom_line(aes(group = cell_line), colour = "grey60") +
  geom_point(aes(colour = treatment), size = 2.5) +
  facet_wrap(~ gene, scales = "free_y", nrow = 2) +
  scale_y_log10() +
  scale_colour_manual(values = TREAT_COLOURS, guide = "none") +
  scale_x_discrete(labels = c(untreated = "ctrl", dexamethasone = "dex")) +
  labs(title = "Known glucocorticoid-responsive genes",
       subtitle = "Normalised counts; lines connect the same cell line",
       x = NULL, y = "Normalised count + 1 (log scale)") +
  theme_report()
save_plot(p, "marker_genes.png", 8, 5.5)

cat("Saved figures to", s$out_dir, "\n")

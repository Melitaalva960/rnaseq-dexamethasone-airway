## =====================================================================
## 05  Compare Phase 1 (published counts) with Phase 2 (your own Salmon run)
##
## Do the two pipelines agree? Compares log2 fold changes for genes tested
## in both, and the overlap of significant genes.
##
## Output: results/figures/phase_comparison.png
##         results/tables/phase_comparison.txt
## =====================================================================
source("R/utils.R")
need("ggplot2")

s <- get_settings(list(
  phase1  = "results/tables/deseq2_all_genes.csv",
  phase2  = "results/phase2/tables/deseq2_all_genes.csv",
  fig     = "results/figures/phase_comparison.png",
  out_txt = "results/tables/phase_comparison.txt"
))
ensure_dir(dirname(s$fig)); ensure_dir(dirname(s$out_txt))

p1 <- read.csv(s$phase1); p2 <- read.csv(s$phase2)
m <- merge(p1[, c("gene_id", "symbol", "log2FC", "padj", "direction")],
           p2[, c("gene_id", "log2FC", "padj", "direction")],
           by = "gene_id", suffixes = c("_p1", "_p2"))

r <- cor(m$log2FC_p1, m$log2FC_p2, use = "complete.obs")
sig1 <- m$gene_id[m$direction_p1 != "not significant"]
sig2 <- m$gene_id[m$direction_p2 != "not significant"]
both <- intersect(sig1, sig2)

lines <- c(
  sprintf("Genes tested in both: %d", nrow(m)),
  sprintf("Pearson correlation of log2 fold changes: %.3f", r),
  sprintf("Significant in Phase 1: %d | Phase 2: %d | both: %d", length(sig1), length(sig2), length(both)),
  sprintf("Of the Phase 1 significant genes, %.1f%% are also significant in Phase 2",
          100 * length(both) / max(1, length(sig1))),
  "", "Marker genes (log2FC Phase 1 vs Phase 2):",
  capture.output(print(m[match(MARKER_GENES, m$symbol), c("symbol", "log2FC_p1", "log2FC_p2")],
                       row.names = FALSE, digits = 3))
)
writeLines(lines, s$out_txt); cat(lines, sep = "\n")

m$class <- ifelse(m$gene_id %in% both, "significant in both",
            ifelse(m$gene_id %in% sig1, "Phase 1 only",
             ifelse(m$gene_id %in% sig2, "Phase 2 only", "neither")))
m$class <- factor(m$class, levels = c("neither", "Phase 1 only", "Phase 2 only", "significant in both"))
p <- ggplot(m[order(m$class), ], aes(log2FC_p1, log2FC_p2, colour = class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 0.7, alpha = 0.7) +
  scale_colour_manual(values = c(neither = "grey80", `Phase 1 only` = "#7A8B99",
                                 `Phase 2 only` = "#FDBA74", `significant in both` = "#C2410C"), name = NULL) +
  labs(title = "Published counts vs my Salmon pipeline",
       subtitle = sprintf("log2 fold change, dexamethasone vs untreated (r = %.2f)", r),
       x = "Phase 1: airway package counts", y = "Phase 2: Salmon + tximport") +
  coord_equal() + theme_report()
ggsave(s$fig, p, width = 7.5, height = 6.5, dpi = 300, bg = "white")
cat("Saved", s$fig, "\n")

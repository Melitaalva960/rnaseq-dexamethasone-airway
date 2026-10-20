## =====================================================================
## 04  Gene Ontology enrichment (clusterProfiler)
##
## Over-representation analysis of GO Biological Processes, run
## separately for up- and down-regulated genes. The background
## ("universe") is every gene that was tested, not the whole genome,
## which avoids inflated results.
##
## Output: results/tables/go_bp_up.csv, go_bp_down.csv
##         results/figures/go_bp_dotplot.png
## =====================================================================
source("R/utils.R")
need(c("clusterProfiler", "org.Hs.eg.db", "ggplot2"))

s <- get_settings(list(
  table      = "results/tables/deseq2_all_genes.csv",
  table_dir  = "results/tables",
  fig_dir    = "results/figures",
  show_terms = 10
))
ensure_dir(s$table_dir); ensure_dir(s$fig_dir)

res <- read.csv(s$table)
universe <- res$gene_id[!is.na(res$padj)]

run_go <- function(genes) {
  ego <- enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                  keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  simplify(ego, cutoff = 0.7)   # merge near-duplicate GO terms
}

out <- list()
for (dir in c("up", "down")) {
  genes <- res$gene_id[res$direction == dir]
  ego <- run_go(genes)
  df <- if (is.null(ego)) data.frame() else as.data.frame(ego)
  write.csv(df, file.path(s$table_dir, sprintf("go_bp_%s.csv", dir)), row.names = FALSE)
  cat(sprintf("%s-regulated: %d genes, %d enriched GO terms\n", dir, length(genes), nrow(df)))
  if (nrow(df) > 0) {
    df <- head(df, s$show_terms)
    df$direction <- if (dir == "up") "Up in dexamethasone" else "Down in dexamethasone"
    out[[dir]] <- df
  }
}

if (length(out) == 0) stop("No enriched GO terms found.")
plot_df <- do.call(rbind, out)
ratio <- function(x) sapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
plot_df$GeneRatio <- ratio(plot_df$GeneRatio)
plot_df$Description <- ifelse(nchar(plot_df$Description) > 55,
                              paste0(substr(plot_df$Description, 1, 52), "..."), plot_df$Description)
plot_df$Description <- factor(plot_df$Description, levels = unique(rev(plot_df$Description)))
plot_df$direction <- factor(plot_df$direction, levels = c("Up in dexamethasone", "Down in dexamethasone"))

p <- ggplot(plot_df, aes(GeneRatio, Description, size = Count, colour = p.adjust)) +
  geom_point() +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
  scale_colour_gradient(low = "#C2410C", high = "#FDBA74", name = "Adjusted p") +
  labs(title = "GO Biological Process enrichment",
       subtitle = sprintf("Top %d terms per direction; background = all tested genes", s$show_terms),
       x = "Gene ratio", y = NULL, size = "Genes") +
  theme_report() +
  theme(strip.text.y = element_text(angle = 0, face = "bold"))
ggsave(file.path(s$fig_dir, "go_bp_dotplot.png"), p, width = 9.5, height = 7.5, dpi = 300, bg = "white")
cat("Saved", file.path(s$fig_dir, "go_bp_dotplot.png"), "\n")

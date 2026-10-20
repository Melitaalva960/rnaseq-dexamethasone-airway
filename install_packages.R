## Install everything the R scripts need (run once):
##   Rscript install_packages.R      or open in RStudio and click "Source"

cran <- c("ggplot2", "ggrepel", "pheatmap", "BiocManager")
bioc <- c("airway", "DESeq2", "apeglm", "tximport", "AnnotationDbi", "org.Hs.eg.db", "clusterProfiler")

missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) install.packages(missing_cran)

missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

ok <- vapply(c(cran, bioc), requireNamespace, logical(1), quietly = TRUE)
print(ok)
if (all(ok)) cat("\nAll packages installed.\n") else cat("\nSome packages failed - see messages above.\n")

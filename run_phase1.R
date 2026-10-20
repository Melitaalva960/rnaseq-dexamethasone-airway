## Phase 1 in one go: counts from the airway package -> DESeq2 -> figures -> GO.
## In RStudio: open the project folder (File > Open Project... or setwd() to it),
## open this file and click "Source".
## From a terminal in the project folder:  Rscript run_phase1.R
##
## Each script uses its default settings (input = airway).
if (!file.exists("R/utils.R")) stop("Run this from the project folder (the one containing R/ and data/).")

for (f in c("R/01_get_counts.R", "R/02_deseq2.R", "R/03_figures.R", "R/04_enrichment.R")) {
  cat("\n==========", f, "==========\n")
  source(f, local = new.env())
}
cat("\nDone. Results are in results/tables and results/figures\n")

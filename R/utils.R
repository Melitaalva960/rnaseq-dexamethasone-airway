## Shared helpers. Every script starts with source("R/utils.R") and is run
## from the repository root, e.g.   Rscript R/02_deseq2.R input=airway

## Settings: each script defines defaults; override any of them as key=value
get_settings <- function(defaults) {
  for (a in commandArgs(trailingOnly = TRUE)) {
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    if (length(kv) < 2 || !kv[1] %in% names(defaults)) {
      stop("Unknown or malformed setting '", a, "'. Valid settings: ",
           paste(names(defaults), collapse = ", "))
    }
    value <- paste(kv[-1], collapse = "=")
    defaults[[kv[1]]] <- if (is.numeric(defaults[[kv[1]]])) as.numeric(value) else value
  }
  defaults
}

ensure_dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

## Stop with a helpful message if a package is missing
need <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing package(s): ", paste(missing, collapse = ", "),
         "\nInstall them by running:  Rscript install_packages.R")
  }
  invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))
}

## Thresholds used throughout
PADJ_MAX <- 0.05   # adjusted p-value (Benjamini-Hochberg)
LFC_MIN  <- 1      # |log2 fold change| >= 1, i.e. at least 2-fold

## Known glucocorticoid-responsive genes: a positive control for the analysis
MARKER_GENES <- c("FKBP5", "TSC22D3", "PER1", "ZBTB16", "KLF15", "DUSP1")

theme_report <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(colour = "grey30"),
                   panel.grid.minor = ggplot2::element_blank())
}

TREAT_COLOURS <- c(untreated = "#7A8B99", dexamethasone = "#C2410C")

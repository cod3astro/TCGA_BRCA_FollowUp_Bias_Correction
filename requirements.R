# ============================================================
# requirements.R
# Run this once before the numbered scripts. Installs everything
# the pipeline needs. Tested with R 4.3.1 and TCGAbiolinks v3.16.1.
# ============================================================

cran_packages <- c("survival", "survminer", "dplyr", "ggplot2")
install.packages(setdiff(cran_packages, rownames(installed.packages())))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"), update = FALSE, ask = FALSE)

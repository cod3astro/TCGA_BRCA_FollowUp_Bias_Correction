# ============================================================
# 01_data_acquisition.R
# Pulls the TCGA-BRCA RNA-seq cohort from the GDC and saves a
# local copy so the rest of the pipeline never has to re-download it.
# Run this first - everything downstream reads data/brca_data.rds.
# ============================================================

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

library(TCGAbiolinks)
library(SummarizedExperiment)

# Query: primary tumor, open-access RNA-seq, STAR-Counts workflow
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  sample.type = "Primary Tumor",
                  access = "open",
                  experimental.strategy = "RNA-Seq")

print(query)
GDCdownload(query)
brca_data <- GDCprepare(query)

saveRDS(brca_data, "data/brca_data.rds")
print("Data successfully saved!")

# Sanity check: confirm it loads back and show what clinical fields we have
brca_data <- readRDS("data/brca_data.rds")
print("Data successfully loaded")
print(dim(brca_data))

clinical_data <- colData(brca_data)
print("Clinical variables available:")
print(colnames(clinical_data)[grep("follow", colnames(clinical_data), ignore.case = TRUE)])

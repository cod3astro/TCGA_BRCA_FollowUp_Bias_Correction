# TCGA-BRCA Survival Analysis: Correcting Follow-Up Time Bias

![R](https://img.shields.io/badge/R-4.3.1-276DC3?logo=r)
![Data](https://img.shields.io/badge/data-TCGA--BRCA-blue)
![Status](https://img.shields.io/badge/status-complete-brightgreen)

Re-analyzing survival outcomes across PAM50 molecular subtypes in the TCGA breast cancer cohort, and correcting a follow-up time bias that was silently making the most aggressive tumor subtype look like it had the best prognosis.

## Table of Contents
- [Background](#background)
- [Data](#data)
- [Method](#method)
- [Key Result](#key-result)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Usage](#usage)
- [Related Writing](#related-writing)
- [Limitations & Next Steps](#limitations--next-steps)
- [Citation](#citation)

## Background

Survival analysis of TCGA's PAM50 molecular subtypes (Luminal A, Luminal B, HER2-enriched, Basal-like, Normal-like) is a standard tool in precision oncology. But multi-institutional datasets like TCGA don't follow every patient for the same length of time, and when follow-up length differs systematically between subgroups, groups with longer observation windows can look like they survive better, even when nothing biological is going on. This project catches that exact bias in the TCGA-BRCA cohort and corrects it.

## Data

| | |
|---|---|
| Source | TCGA-BRCA (Breast Invasive Carcinoma), via `TCGAbiolinks` |
| Query | Transcriptome Profiling · Gene Expression Quantification · STAR-Counts · Primary Tumor · open access |
| Cohort size | 1,111 patients queried -> 1,107 retained (4 excluded for missing survival data) |

**PAM50 subtype breakdown**

| Subtype | n |
|---|---|
| Luminal A | 567 |
| Luminal B | 209 |
| HER2-enriched | 82 |
| Basal-like | 197 |
| Normal-like | 40 |

## Method

1. Build a survival object from `days_to_death` (deceased) / `days_to_last_follow_up` (censored)
2. Fit Kaplan-Meier curves (`survfit()`) and compare subtypes with log-rank tests (`survdiff()`)
3. Check maximum follow-up length per subtype -> find large disparities (4,267 to 8,605 days)
4. Truncate every patient to the shortest common observation window (4,267 days) and re-run the comparison

## Key Result

| Subtype | Median OS, uncorrected (days) | Median OS, corrected (days) | Rank shift |
|---|---|---|---|
| Basal-like | 7,455 | 3,472 | 1st -> 4th |
| HER2-enriched | 6,456 | 3,063 | 2nd -> 5th |
| Normal-like | 4,267 | 4,267 | unchanged (1st) |
| Luminal A | 3,945 | 3,462 | roughly unchanged |
| Luminal B | 3,941 | 3,941 | roughly unchanged |

Before correction, Basal-like, the clinically most aggressive subtype, appeared to have the *best* survival, purely because those patients had the longest follow-up window in the dataset. After truncating every subtype to the same 4,267-day window, the ranking flipped to match known breast cancer biology, with HER2-enriched correctly showing the poorest survival. The corrected five-way comparison was borderline significant (chi-sq = 8.1, df = 4, p = 0.09) consistent with small subtype sample sizes, not an absence of real difference.

## Repository Structure

```
├── scripts/
│   ├── 01_data_acquisition.R
│   ├── 02_survival_analysis.R
│   ├── 03_bias_correction.R
│   └── 04_figures.R
├── data/
│   └── gdc_manifest.txt     # moved out of the repo root
├── results/
├── Figures/
├── requirements.R           # or renv.lock
└── LICENSE
```

## Requirements

```r
install.packages(c("survival", "survminer", "dplyr", "ggplot2"))
# Bioconductor packages:
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"))
```
Tested with R 4.3.1 and TCGAbiolinks v3.16.1.

## Usage

```bash
git clone https://github.com/cod3astro/project-oncoSurvive.git
cd project-oncoSurvive
```
```r
source("scripts/analysis.R")
```
Note: `GDCdownload()` will pull the full TCGA-BRCA RNA-seq cohort (several GB) on first run.

## Related Writing

Full write-up on Medium: [Methodological Considerations in TCGA Breast Cancer Survival Analysis](https://medium.com/@cod3astro/methodological-considerations-in-tcga-breast-cancer-survival-analysis-6f8bf9b3dce7)

## Limitations & Next Steps

- The truncation approach trades statistical power for a fair comparison — restricted mean survival time or time-dependent covariate models could recover some of that power without reintroducing the bias.
- Small subtype sample sizes (HER2-enriched n=82, Normal-like n=40) limit how confidently the corrected ranking can be generalized.

## Citation

If you use this analysis or approach, please cite the write-up above. TCGA data citation: Cancer Genome Atlas Network (2012), *Nature* 490:61-70.

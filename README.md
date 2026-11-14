# TCGA-BRCA Survival Analysis: Follow-Up Time and Median Survival

![R](https://img.shields.io/badge/R-4.5-blue) ![Status](https://img.shields.io/badge/status-corrected%20reanalysis-green)

An R reanalysis of survival across PAM50 breast cancer subtypes in the TCGA-BRCA cohort. The project asks one question: when patients are followed for very different lengths of time, which survival numbers can be trusted, and which depend on a handful of patients at the tail of the curve?

## Table of Contents

- [Background](#background)
- [Data](#data)
- [Method](#method)
- [Key Results](#key-results)
- [What Was Corrected](#what-was-corrected)
- [Limitations](#limitations)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Usage](#usage)
- [Related Writing](#related-writing)
- [Citation](#citation)

## Background

Survival of the PAM50 subtypes (Luminal A, Luminal B, HER2-enriched, Basal-like, Normal-like) is a standard comparison in breast cancer research. In TCGA, patients were not followed for the same length of time, and the follow-up of each subtype ends at a different day.

A Kaplan-Meier curve handles different follow-up lengths correctly in principle. The risk is at the end of each curve, where very few patients are still being followed. A median survival time is read from the point where the curve crosses 50%. When few patients died, that crossing can happen only at the tail, where one death among a few remaining patients moves the curve a long way. A median read there describes those few patients, not the subtype.

## Data

| | |
|---|---|
| Source | TCGA-BRCA (Breast Invasive Carcinoma), through `TCGAbiolinks` |
| Query | Transcriptome Profiling, Gene Expression Quantification, STAR-Counts, Primary Tumor, open access |
| Samples queried | 1,111 primary-tumour samples |
| After cleaning | 1,107 samples had a usable survival time. Five of these were repeat samples from a patient already counted, so only the first sample per patient was kept: **1,102 unique patients** |
| With a PAM50 call | 1,090 patients (Basal-like 195, HER2-enriched 82, Luminal A 564, Luminal B 209, Normal-like 40) |

## Method

1. Build the survival time: days to death for patients who died, days to last follow-up for the rest (censored).
2. Fit Kaplan-Meier curves for the whole cohort and by subtype. Count how many patients are still at risk at late time points.
3. Find the shortest follow-up among the subtypes. Use it as a shared window for all subtypes.
4. Apply the window by censoring, never by deleting. A patient followed past the window stays in the data and counts as alive at the window edge. Deleting those patients would remove the longest survivors and make early survival look worse.
5. Compare the subtypes inside the window with: the share of patients alive at 5 years, restricted mean survival time (average days alive inside the window), a log-rank test, and a Cox model with Luminal A as the reference. Check the Cox model's proportional-hazards assumption.
6. Repeat the comparison without the Normal-like subtype, because the smallest subtype sets the window for all the others.

Medians are still reported for comparison, but they are not used for conclusions. "Not reached" means the curve never fell below 50% inside the window.

## Key Results

Every number below is written by `scripts/03_bias_correction.R` to `results/results_summary.md`.

Shared window: **4,267 days**, set by the Normal-like subtype (shortest follow-up).
Log-rank test across subtypes inside the window: **p = 0.0345**.
Same test without Normal-like (window 4,894 days): **p = 0.0198**.

### Patients, deaths and follow-up (days)

| subtype | patients | deaths_total | deaths_in_window | longest_followup |
|---|---|---|---|---|
| Basal | 195 | 28 | 27 | 8605 |
| Her2 | 82 | 17 | 16 | 6456 |
| LumA | 564 | 66 | 64 | 8556 |
| LumB | 209 | 33 | 33 | 4894 |
| Normal | 40 | 7 | 7 | 4267 |

### Share alive at 5 years inside the window

| subtype | at_risk | alive_5y | CI_95 |
|---|---|---|---|
| Basal | 47 | 0.791 | 0.713 to 0.877 |
| Her2 | 16 | 0.718 | 0.595 to 0.865 |
| LumA | 134 | 0.867 | 0.825 to 0.912 |
| LumB | 41 | 0.769 | 0.681 to 0.869 |
| Normal | 9 | 0.674 | 0.477 to 0.953 |

### Average days alive inside the window (restricted mean survival time)

| subtype | rmst_days | std_error |
|---|---|---|
| Basal | 3375 | 152 |
| Her2 | 3012 | 269 |
| LumA | 3397 | 91 |
| LumB | 3000 | 163 |
| Normal | 3337 | 327 |

### Cox model: hazard ratio against Luminal A (inside the window)

| subtype | HR | CI_95 | p_value |
|---|---|---|---|
| Basal | 1.21 | 0.77 to 1.89 | 0.41 |
| Her2 | 2.11 | 1.22 to 3.65 | 0.0077 |
| LumB | 1.64 | 1.08 to 2.50 | 0.021 |
| Normal | 1.46 | 0.67 to 3.19 | 0.34 |

Proportional-hazards check: global p = 0.017.

### Patients still followed at late time points (before any correction)

| subtype | day 3000 | day 4000 | day 5000 | day 6000 | day 7000 |
|---|---|---|---|---|---|
| Basal | 22 | 7 | 4 | 4 | 4 |
| Her2 | 5 | 2 | 1 | 1 | 0 |
| LumA | 48 | 18 | 11 | 6 | 4 |
| LumB | 13 | 2 | 0 | 0 | 0 |
| Normal | 4 | 1 | 0 | 0 | 0 |

### Median survival in days (not reached = curve never fell below 50%)

| subtype | median_uncorrected | median_in_window |
|---|---|---|
| Basal | 7455 | not reached |
| Her2 | 6456 | not reached |
| LumA | 3945 | 3945 |
| LumB | 3941 | 3941 |
| Normal | 4267 | 4267 |

### How to read these results

1. **The uncorrected medians rest on a few patients.** The Basal-like median of 7,455 days is the day of the last death, when only four Basal-like patients were still followed. HER2-enriched (6,456 days) and Normal-like (4,267 days) are the same: each median is the last event of a nearly empty tail. They do not show that Basal-like patients survive longer.
2. **The subtypes differ inside the shared window** (log-rank p = 0.0345, and 0.0198 without Normal-like).
3. **HER2-enriched and Luminal B carry a higher risk of death than Luminal A** (hazard ratios 2.11 and 1.64, confidence intervals above 1). Their average days alive inside the window are about 400 days lower than Luminal A's.
4. **Basal-like is not clearly different from Luminal A** in the Cox model (HR 1.21, interval 0.77 to 1.89).
5. **Read the hazard ratios with caution.** The proportional-hazards check is significant (p = 0.017), so one hazard ratio per subtype does not describe the whole follow-up. The estimates hint that Basal-like does worse than Luminal A early (alive at 5 years: 0.791 vs 0.867) and similarly or better later (alive at 10 years: 0.633 vs 0.589, but only 8 and 24 patients are still followed, and the intervals overlap widely). Restricted mean survival time and the 5-year figures are the safer summaries.

## What Was Corrected

The first version of this project contained three errors, found and fixed in November 2026:

1. **Deleted patients.** The window was applied by removing patients followed longer than the window. They should have been censored at the window edge. Deleting them removed long survivors and lowered early survival estimates.
2. **Hand-typed results.** The before-and-after table was typed into the code, and three of its "corrected medians" (Basal-like 3,472, HER2-enriched 3,063, Luminal A 3,462) were the lower confidence limits of the uncorrected fit, not medians. All tables are now written by the code to `results/results_summary.md`.
3. **Overstated conclusions.** The first write-up described a reversed subtype ranking, a non-significant comparison (p = 0.09) and validated survival curves. The corrected analysis does not support that wording. The Normal-like subtype was also described as normal tissue, but it is a tumour label.

Five repeat samples from patients already counted were also removed. The earlier counts included them.

## Limitations

- About 150 deaths in total, and at most 66 in any one subtype. The 5-year estimates have wide confidence intervals, especially for HER2-enriched (16 patients at risk) and Normal-like (9 patients at risk).
- Estimates at 10 years rest on very few patients (2 for HER2-enriched, 1 for Normal-like) and are not used as findings.
- The proportional-hazards assumption of the Cox model is questionable (check p = 0.017). Hazard ratios are an average over follow-up, not a fixed effect.
- No adjustment for age, stage, treatment or other factors. The subtype comparison is descriptive, not causal.
- The shared window is set by the smallest subtype (Normal-like, 40 patients). This is why the comparison is repeated without it. The answer does not change.
- Several statistics were computed. The log-rank p-values of 0.03 and 0.02 are modest evidence, not proof.
- Time-dependent or stratified models are a possible extension.

## Repository Structure

```
├── scripts/
│   ├── 01_data_acquisition.R    # download and save the TCGA-BRCA cohort
│   ├── 02_survival_analysis.R   # survival object, Kaplan-Meier fits, at-risk check
│   ├── 03_bias_correction.R     # shared window, censoring, tests, results_summary.md
│   ├── 04_figures.R             # every figure, session info
│   └── run_all.R                # runs the four scripts in order
├── data/                        # downloaded data (not tracked, several GB)
├── results/                     # saved objects, results_summary.md, session_info.txt
└── Figures/                     # all plots, with risk tables where relevant
```

## Requirements

```r
install.packages(c("survival", "survminer", "dplyr", "ggplot2", "ggpubr"))
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"))
```

Tested with R 4.5.1, survival 3.8.3, survminer 0.5.1 and TCGAbiolinks 2.36.0. Full session details are saved in `results/session_info.txt`.

## Usage

Set the working directory to the project folder, not to `scripts/`, then run:

```r
source("scripts/run_all.R")
```

The first run downloads about 4.7 GB of RNA-seq data and takes a while. Later runs skip the download if `data/brca_data.rds` exists. To run a single script, run the ones before it first, because each script reads files saved by the earlier ones.

## Related Writing

The earlier Medium write-up, *Methodological Considerations in TCGA Breast Cancer Survival Analysis*, describes the first version of the method. See [What Was Corrected](#what-was-corrected) before relying on it.

## Citation

If you use this analysis or approach, please cite this repository. TCGA data citation: Cancer Genome Atlas Network (2012), Nature 490:61-70.

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
| Cleaning | Samples without a usable survival time are dropped. Some patients appear under more than one sample, so only the first sample per patient is kept. Samples without a PAM50 call are dropped for the subtype analysis. |

Final patient counts per subtype are in the results block below, produced by the code.

## Method

1. Build the survival time: days to death for patients who died, days to last follow-up for the rest (censored).
2. Fit Kaplan-Meier curves for the whole cohort and by subtype. Count how many patients are still at risk at late time points.
3. Find the shortest follow-up among the subtypes. Use it as a shared window for all subtypes.
4. Apply the window by censoring, never by deleting. A patient followed past the window stays in the data and counts as alive at the window edge. Deleting those patients would remove the longest survivors and make early survival look worse.
5. Compare the subtypes inside the window with: the share of patients alive at 5 years, restricted mean survival time (average days alive inside the window), a log-rank test, and a Cox model with Luminal A as the reference.
6. Repeat the comparison without the Normal-like subtype, because the smallest subtype may be the one that sets the window for all the others.

Medians are still reported for comparison, but they are not used for conclusions. A median marked "not reached" means the curve never fell below 50% inside the window.

## Key Results

<!-- Paste the contents of results/results_summary.md here after running the pipeline, then delete this line. -->

**How to read the results**

- Compare subtypes using the 5-year survival table and the Cox hazard ratios, not the medians. A hazard ratio above 1 means a higher risk of death than Luminal A.
- The late-time at-risk table shows why uncorrected medians are fragile: past about 4,000 days only a few patients remain in some subtypes.
- Check the three summary sentences below against the tables before keeping them:
  1. The subtypes differ in survival inside the shared window (log-rank).
  2. HER2-enriched and Luminal B have a higher risk than Luminal A.
  3. Basal-like is not clearly different from Luminal A, and the high uncorrected Basal-like median came from the tail, not from better survival.

## What Was Corrected

The first version of this project contained three errors, found and fixed in November 2026:

1. **Deleted patients.** The window was applied by removing patients followed longer than the window. They should have been censored at the window edge. Deleting them removed long survivors and lowered early survival estimates.
2. **Hand-typed results.** The before-and-after table was typed into the code, and two of its "corrected medians" matched the lower confidence limits of the uncorrected fit, not medians. All tables are now written by the code to `results/results_summary.md`.
3. **Overstated conclusions.** The first write-up described a reversed subtype ranking and called the survival curves validated. The corrected analysis does not support that wording. The Normal-like subtype was also described as normal tissue, but it is a tumour label.

Repeat samples from the same patient were also removed. The earlier counts included them.

## Limitations

- About 150 deaths in total, and fewer than 70 in any one subtype. The 5-year estimates have wide confidence intervals, especially for HER2-enriched and Normal-like.
- Estimates at 10 years rest on very few patients and are not reported as findings.
- No adjustment for age, stage, treatment or other factors. The subtype comparison is descriptive, not causal.
- The shared window is set by the subtype with the shortest follow-up, which here is the smallest group. This is why the comparison is repeated without it.
- The Cox model assumes risks stay in proportion over time. The script prints a check, and the hazard ratios should be read with that check in mind.
- Restricted mean survival time and time-dependent methods are possible extensions.

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

Tested with R 4.5.1. Exact package versions from the last run are saved in `results/session_info.txt`.

## Usage

Set the working directory to the project folder, not to `scripts/`, then run:

```r
source("scripts/run_all.R")
```

The first run downloads several GB of RNA-seq data and takes a while. Later runs skip the download if `data/brca_data.rds` exists. To run a single script, run the ones before it first, because each script reads files saved by the earlier ones.

## Related Writing

The earlier Medium write-up, *Methodological Considerations in TCGA Breast Cancer Survival Analysis*, describes the first version of the method. See [What Was Corrected](#what-was-corrected) before relying on it.

## Citation

If you use this analysis or approach, please cite this repository. TCGA data citation: Cancer Genome Atlas Network (2012), Nature 490:61-70.

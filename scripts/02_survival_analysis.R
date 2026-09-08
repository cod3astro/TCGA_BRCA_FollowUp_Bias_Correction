# ============================================================
# 02_survival_analysis.R
# Builds the survival object from the raw clinical data and fits
# baseline Kaplan-Meier curves - overall cohort, then by PAM50 subtype.
# This is the "naive" analysis, before we check it for bias in 03.
# Needs data/brca_data.rds from 01_data_acquisition.R.
# ============================================================

library(survival)
library(survminer)
library(dplyr)
library(SummarizedExperiment)

brca_data <- readRDS("data/brca_data.rds")
clinical_data <- colData(brca_data)

# --- Build the survival dataframe ---
survival_df <- data.frame(
  patient_id = clinical_data$bcr_patient_barcode,
  vital_status = clinical_data$vital_status,
  days_to_death = clinical_data$days_to_death,
  days_to_last_followup = clinical_data$days_to_last_follow_up
)

print(table(survival_df$vital_status))

# time-to-event: days to death if deceased, else days to last follow-up
survival_df$time <- ifelse(survival_df$vital_status == "Dead",
                            survival_df$days_to_death,
                            survival_df$days_to_last_followup)
survival_df$event <- ifelse(survival_df$vital_status == "Dead", 1, 0)

# Drop patients with no usable time value
survival_df <- survival_df[!is.na(survival_df$time), ]
print(paste("Patients for survival analysis:", nrow(survival_df)))

# --- Overall Kaplan-Meier fit ---
surv_object <- Surv(time = survival_df$time, event = survival_df$event)
km_fit <- survfit(surv_object ~ 1, data = survival_df)

summary_stats <- summary(km_fit)
median_survival <- summary_stats$table["median"]
print("Survival summary statistics:")
print(paste("Median survival time:", round(median_survival, 2), "days"))
print(paste("Number of patients:", km_fit$n))
print(paste("Number of events (death):", sum(survival_df$event)))

time_points <- c(365, 739, 1825)
surv_at_times <- summary(km_fit, time = time_points)
print("Survival probabilities at key time points:")
print(surv_at_times$surv)

max_observed_time <- max(summary(km_fit)$time)
print(paste("Maximum observed follow-up time:", max_observed_time, "days"))

# --- Add PAM50 subtype and fit per-subtype curves ---
print("Available subtype columns:")
print(colnames(clinical_data)[grep("subtype|pam50|brca_subtype", colnames(clinical_data), ignore.case = TRUE)])

survival_df$subtype <- clinical_data$paper_BRCA_Subtype_PAM50[match(
  survival_df$patient_id, clinical_data$bcr_patient_barcode)]

# Drop patients with no subtype call
subtype_survival <- survival_df[!is.na(survival_df$subtype), ]
print("Patients with subtype information:")
print(table(subtype_survival$subtype))

subtype_surv_object <- Surv(time = subtype_survival$time, event = subtype_survival$event)
km_subtype_fit <- survfit(subtype_surv_object ~ subtype, data = subtype_survival)

print("Median survival by subtype (uncorrected):")
print(surv_median(km_subtype_fit))
print("Number of deaths by subtype:")
print(table(subtype_survival$subtype, subtype_survival$event))

# --- Save everything 03 and 04 need, so they don't have to refit ---
saveRDS(survival_df, "results/survival_df.rds")
saveRDS(subtype_survival, "results/subtype_survival.rds")
saveRDS(km_fit, "results/km_overall.rds")
saveRDS(km_subtype_fit, "results/km_subtype.rds")
print("Saved survival_df, subtype_survival, km_overall, km_subtype to results/")

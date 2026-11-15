# ============================================================
# 03_bias_correction.R
# Checks whether follow-up length differs between PAM50 subtypes,
# then compares the subtypes inside one shared time window so each
# one is measured over the same stretch of time.
#
# The rule that matters most in this file: a patient who is still
# alive past the window edge is CENSORED at the edge, never deleted.
# Deleting them removes the longest survivors and makes survival look
# worse than it really was.
#
# Needs results/subtype_survival.rds from 02_survival_analysis.R.
# Writes the .rds files that 04_figures.R reads, and results/results_summary.md
# with every number the README quotes, so none are copied by hand.
# ============================================================

library(survival)
library(survminer)
library(dplyr)

subtype_survival <- readRDS("results/subtype_survival.rds")

# ------------------------------------------------------------
# cut_at_window(): puts every patient inside the same time window.
# Takes the survival table and a cut-off in days. Gives back a table
# with the SAME patients, where anyone observed past the cut-off now
# counts as "still alive at the cut-off" (event = 0, time = cut-off).
# Why it exists: a fair comparison needs a shared window, but the first
# version of this project deleted everyone followed longer than the
# window, which removed exactly the patients who survived the longest.
# ------------------------------------------------------------
cut_at_window <- function(df, cutoff) {
  df %>%
    mutate(
      # event has to be changed first, because it needs the original time
      event = ifelse(time > cutoff, 0, event),
      time  = pmin(time, cutoff)
    )
}

# --- Step 1: how different is follow-up between subtypes? ---
followup_check <- subtype_survival %>%
  group_by(subtype) %>%
  summarize(
    n_patients      = n(),
    deaths          = sum(event),
    max_followup    = max(time),
    median_followup = median(time),
    n_beyond_4000   = sum(time > 4000)  # how thin the late tail is
  )
print("Follow-up length and number of deaths by subtype (before any correction):")
print(followup_check)

# --- Step 2: choose the shared window and apply it ---
common_followup <- min(followup_check$max_followup)
window_setter   <- followup_check$subtype[which.min(followup_check$max_followup)]
print(paste("Shared window:", common_followup, "days. It is set by the subtype:", window_setter))

truncated_survival <- cut_at_window(subtype_survival, common_followup)

# Hard stop if a patient went missing. The old version silently lost 23.
stopifnot(nrow(truncated_survival) == nrow(subtype_survival))
print("Patients per subtype after applying the window (must equal the full counts):")
print(table(truncated_survival$subtype))

# --- Step 3: refit Kaplan-Meier inside the window ---
km_trunc_fit <- survfit(Surv(time, event) ~ subtype, data = truncated_survival)

# --- Step 4: measures that do not depend on one fragile median ---
# Medians only exist where a curve drops below 50%, which with few deaths
# happens late, where few patients are left. So report three other things too.

# 4a. Share of patients alive at 5 and 10 years, with confidence limits
alive_at <- summary(km_trunc_fit, times = c(1825, 3650), extend = TRUE)
survival_at_times <- data.frame(
  subtype  = sub("subtype=", "", as.character(alive_at$strata)),
  days     = alive_at$time,
  n_risk   = alive_at$n.risk,
  survival = round(alive_at$surv, 3),
  lower    = round(alive_at$lower, 3),
  upper    = round(alive_at$upper, 3)
)
print("Share of patients alive at 5 years (1825 d) and 10 years (3650 d):")
print(survival_at_times)

# 4b. Restricted mean survival time: average days lived inside the window
rmst_table <- summary(km_trunc_fit, rmean = common_followup)$table
print("Restricted mean survival time (average days alive inside the window):")
print(rmst_table[, c("records", "events", "rmean", "se(rmean)")])

# 4c. Medians before and after, computed from the fits (NA = curve never fell below 50%)
original_fit <- survfit(Surv(time, event) ~ subtype, data = subtype_survival)
orig_med  <- surv_median(original_fit)
trunc_med <- surv_median(km_trunc_fit)
comparison <- data.frame(
  subtype         = sub("subtype=", "", orig_med$strata),
  median_original = orig_med$median,
  median_window   = trunc_med$median
)
print("Median survival before and after the shared window (NA means never reached 50%):")
print(comparison)

# --- Step 5: are the subtypes actually different? ---
subtype_test <- survdiff(Surv(time, event) ~ subtype, data = truncated_survival)
p_logrank <- 1 - pchisq(subtype_test$chisq, df = length(subtype_test$n) - 1)
print("Log-rank test across all subtypes, inside the shared window:")
print(subtype_test)
print(paste("Log-rank p-value:", signif(p_logrank, 3)))

# Same test on the uncorrected data, to show whether the shared window created the
# difference between subtypes or only removed the unstable tail.
uncorrected_test <- survdiff(Surv(time, event) ~ subtype, data = subtype_survival)
p_uncorrected <- 1 - pchisq(uncorrected_test$chisq, df = length(uncorrected_test$n) - 1)
print(paste("Log-rank p-value before any window (uncorrected):", signif(p_uncorrected, 3)))

basal_vs_her2 <- survdiff(Surv(time, event) ~ subtype,
                          data = truncated_survival[truncated_survival$subtype %in% c("Basal", "Her2"), ])
print("Basal vs HER2 only, inside the shared window:")
print(basal_vs_her2)

# Cox model gives one hazard ratio per subtype against Luminal A.
# Works on a copy so the factor order in truncated_survival stays untouched for 04.
cox_data <- truncated_survival
cox_data$subtype <- relevel(factor(cox_data$subtype), ref = "LumA")
cox_fit <- coxph(Surv(time, event) ~ subtype, data = cox_data)
print("Cox model, hazard ratio of each subtype vs Luminal A:")
print(summary(cox_fit))

# The Cox model assumes each subtype's risk stays in proportion over time.
# A small p-value here means the hazard ratios above should be read with caution.
print("Proportional hazards check (small p-value = the assumption may not hold):")
print(cox.zph(cox_fit))

# --- Step 6: does the answer change without the tiniest subtype? ---
# The shared window is set by whichever subtype has the shortest follow-up.
# If that is the very small Normal-like group, it shrinks the window for
# every other subtype too. Repeat the comparison without it.
others <- subtype_survival %>% filter(subtype != "Normal")
window_without_normal <- min(tapply(others$time, others$subtype, max))
print(paste("Window without Normal-like:", window_without_normal, "days"))

sens_data <- cut_at_window(others, window_without_normal)
stopifnot(nrow(sens_data) == nrow(others))
sens_fit  <- survfit(Surv(time, event) ~ subtype, data = sens_data)
sens_test <- survdiff(Surv(time, event) ~ subtype, data = sens_data)
sens_p    <- 1 - pchisq(sens_test$chisq, df = length(sens_test$n) - 1)
print("Median survival without Normal-like (NA means never reached 50%):")
print(surv_median(sens_fit))
print(paste("Log-rank p-value without Normal-like:", signif(sens_p, 3)))

# --- Step 7: what the results say, written from the numbers above ---
print(paste("Corrected log-rank p-value:", signif(p_logrank, 3)))
if (p_logrank < 0.05) {
  print("Below 0.05: the subtypes differ in survival inside the shared window. Read the Cox hazard ratios to see which ones.")
} else {
  print("0.05 or above: the data do not show a clear survival difference between subtypes inside the shared window.")
}
print("Note: Normal-like is a tumour label that resembles normal breast tissue. It is not healthy tissue.")

# --- Step 8: write the numbers the README quotes, so none are copied by hand ---

# ------------------------------------------------------------
# fmt(): fixed-decimal text for table cells. Takes numbers and a digit
# count, gives back strings, so the tables never show 0.7850000001.
# ------------------------------------------------------------
fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

# ------------------------------------------------------------
# md_table(): turns a data frame into a markdown pipe table.
# Takes a data frame, gives back one string per line, ready for writeLines().
# Why it exists: the README tables should come straight from the objects
# above, because the first version of this project had hand-typed numbers.
# ------------------------------------------------------------
md_table <- function(df) {
  df[] <- lapply(df, as.character)
  header  <- paste0("| ", paste(names(df), collapse = " | "), " |")
  divider <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  rows    <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(header, divider, unname(rows))
}

deaths_in_window <- truncated_survival %>%
  group_by(subtype) %>%
  summarize(deaths_in_window = sum(event))

cohort_tbl <- data.frame(
  subtype          = followup_check$subtype,
  patients         = followup_check$n_patients,
  deaths_total     = followup_check$deaths,
  deaths_in_window = deaths_in_window$deaths_in_window,
  longest_followup = followup_check$max_followup
)

five_year <- survival_at_times[survival_at_times$days == 1825, ]
five_year_tbl <- data.frame(
  subtype   = five_year$subtype,
  at_risk   = five_year$n_risk,
  alive_5y  = fmt(five_year$survival),
  CI_95     = paste0(fmt(five_year$lower), " to ", fmt(five_year$upper))
)

cox_sum <- summary(cox_fit)
cox_tbl <- data.frame(
  subtype = sub("subtype", "", rownames(cox_sum$conf.int)),
  HR      = fmt(cox_sum$conf.int[, "exp(coef)"], 2),
  CI_95   = paste0(fmt(cox_sum$conf.int[, "lower .95"], 2), " to ",
                   fmt(cox_sum$conf.int[, "upper .95"], 2)),
  p_value = signif(cox_sum$coefficients[, "Pr(>|z|)"], 2),
  row.names = NULL
)

# How many patients are still being followed late on, before any correction.
# This is the evidence for why a median there rests on very few people.
late_tail <- summary(original_fit, times = c(3000, 4000, 5000, 6000, 7000), extend = TRUE)
late_long <- data.frame(
  subtype = sub("subtype=", "", as.character(late_tail$strata)),
  days    = late_tail$time,
  at_risk = late_tail$n.risk
)
late_wide <- reshape(late_long, idvar = "subtype", timevar = "days", direction = "wide")
names(late_wide) <- sub("at_risk.", "day ", names(late_wide), fixed = TRUE)
rownames(late_wide) <- NULL

rmst_md_tbl <- data.frame(
  subtype   = sub("subtype=", "", rownames(rmst_table)),
  rmst_days = round(rmst_table[, "rmean"]),
  std_error = round(rmst_table[, "se(rmean)"]),
  row.names = NULL
)

median_tbl <- data.frame(
  subtype            = comparison$subtype,
  median_uncorrected = ifelse(is.na(comparison$median_original), "not reached",
                              as.character(comparison$median_original)),
  median_in_window   = ifelse(is.na(comparison$median_window), "not reached",
                              as.character(comparison$median_window))
)

summary_lines <- c(
  "## Results written by 03_bias_correction.R",
  "",
  paste0("Unique patients with a PAM50 subtype: ", nrow(subtype_survival)),
  paste0("Shared window: ", common_followup, " days (set by the subtype with the shortest follow-up: ", window_setter, ")"),
  paste0("Log-rank test across subtypes before any window (uncorrected): p = ", signif(p_uncorrected, 3)),
  paste0("Log-rank test across subtypes inside the window: p = ", signif(p_logrank, 3)),
  paste0("Same test without Normal-like (window ", window_without_normal, " days): p = ", signif(sens_p, 3)),
  "",
  "### Patients, deaths and follow-up (days)", "", md_table(cohort_tbl), "",
  "### Share alive at 5 years inside the window", "", md_table(five_year_tbl), "",
  "### Average days alive inside the window (restricted mean survival time)", "", md_table(rmst_md_tbl), "",
  "### Cox model: hazard ratio against Luminal A (inside the window)", "", md_table(cox_tbl), "",
  "### Patients still followed at late time points (before any correction)", "", md_table(late_wide), "",
  "### Median survival in days (not reached = curve never fell below 50%)", "", md_table(median_tbl), ""
)
writeLines(summary_lines, "results/results_summary.md")
print("Wrote results/results_summary.md. Paste its tables into the README under Key Results.")

# --- Save what 04_figures.R needs ---
saveRDS(truncated_survival, "results/truncated_survival.rds")
saveRDS(km_trunc_fit,       "results/km_truncated.rds")
saveRDS(sens_data,          "results/sensitivity_survival.rds")
saveRDS(sens_fit,           "results/km_sensitivity.rds")
saveRDS(survival_at_times,  "results/survival_at_times.rds")
saveRDS(comparison,         "results/comparison_table.rds")
print("Saved truncated_survival, km_truncated, sensitivity data and fit, survival_at_times, comparison_table to results/")

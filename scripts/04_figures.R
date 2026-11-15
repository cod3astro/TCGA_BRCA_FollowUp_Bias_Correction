# ============================================================
# 04_figures.R
# Regenerates every figure in Figures/ from the saved results,
# so the plots are reproducible from a script instead of a manual
# RStudio export. Needs the .rds files written by 02 and 03.
# The subtype figures are saved with their risk table underneath,
# because the risk table shows how few patients are left late on.
# ============================================================

library(survival)
library(survminer)
library(ggplot2)

dir.create("Figures", showWarnings = FALSE)

# ggpubr (used by survminer) still calls a ggplot2 feature that was renamed,
# which prints a harmless deprecation warning. The plots are not affected.
options(lifecycle_verbosity = "quiet")

survival_df        <- readRDS("results/survival_df.rds")
subtype_survival   <- readRDS("results/subtype_survival.rds")
truncated_survival <- readRDS("results/truncated_survival.rds")
sens_data          <- readRDS("results/sensitivity_survival.rds")
km_fit             <- readRDS("results/km_overall.rds")
km_subtype_fit     <- readRDS("results/km_subtype.rds")
km_trunc_fit       <- readRDS("results/km_truncated.rds")
km_sens_fit        <- readRDS("results/km_sensitivity.rds")

# --- 1. Basic overall survival curve ---
overall_plot <- ggsurvplot(km_fit,
                            data = survival_df,
                            risk.table = TRUE,
                            title = "Overall Survival - TCGA BRCA Patients",
                            xlab = "Time (Days)",
                            ylab = "Survival Probability",
                            legend = "none")
suppressMessages(ggsave("Figures/overall-survival.png", plot = overall_plot$plot, width = 8, height = 6, dpi = 150))

# --- 2. Enhanced overall survival curve (CI, p-value, median line) ---
enhanced_plot <- suppressMessages(ggsurvplot(km_fit,
                             data = survival_df,
                             risk.table = TRUE,
                             conf.int = TRUE,
                             title = "Overall Survival Analysis - TCGA BRCA Cohort",
                             xlab = "Time (Days)",
                             ylab = "Survival Probability",
                             legend = "none",
                             ggtheme = theme_minimal(),
                             risk.table.height = 0.25,
                             surv.median.line = "hv",
                             tables.theme = theme_cleantable()))
enhanced_plot$plot <- enhanced_plot$plot +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"))
suppressMessages(ggsave("Figures/enhanced_survival_plot.png", plot = enhanced_plot$plot, width = 8, height = 6, dpi = 150))

# --- 3. Observed-only survival curve (no extrapolation past last event) ---
max_observed_time <- max(summary(km_fit)$time)
observed_plot <- suppressMessages(ggsurvplot(km_fit,
                             data = survival_df,
                             risk.table = TRUE,
                             conf.int = TRUE,
                             xlim = c(0, max_observed_time),
                             title = "Observed Survival (No Extrapolation) - TCGA BRCA",
                             xlab = "Time (Days)",
                             ylab = "Survival Probability",
                             legend = "none"))
suppressMessages(ggsave("Figures/observed-survival.png", plot = observed_plot$plot, width = 8, height = 6, dpi = 150))

# --- 4. Survival by subtype, uncorrected ---
# Saved with png()/print() instead of ggsave() so the risk table is kept.
subtype_plot <- ggsurvplot(km_subtype_fit,
                            data = subtype_survival,
                            risk.table = TRUE,
                            conf.int = FALSE,
                            pval = TRUE,
                            title = "Survival by Breast Cancer Subtypes - TCGA BRCA",
                            xlab = "Time (Days)",
                            ylab = "Survival Probability",
                            legend.title = "PAM50 Subtype",
                            legend.labs = sub("subtype=", "", names(km_subtype_fit$strata)),
                            risk.table.height = 0.35,
                            risk.table.fontsize = 3.5,
                            break.time.by = 1000)
png("Figures/survival-by-brca-subtypes.png", width = 1200, height = 1300, res = 150)
suppressMessages(print(subtype_plot))
dev.off()

# --- 5. Follow-up time and event boxplot by subtype (the bias itself, visualized) ---
followup_box <- ggplot(subtype_survival, aes(x = subtype, y = time, color = as.factor(event))) +
  geom_boxplot() +
  labs(title = "Follow-up Time and Events by Subtype",
       x = "Subtype",
       y = "Time (Days)",
       color = "Event (0=Alive, 1=Dead)") +
  theme_minimal()
suppressMessages(ggsave("Figures/follow-up-time-and-event-by-subtype-boxplot.png", plot = followup_box, width = 8, height = 6, dpi = 150))

# --- 6. Survival by subtype, inside the shared window (patients censored, not removed) ---
trunc_plot <- ggsurvplot(km_trunc_fit,
                          data = truncated_survival,
                          risk.table = TRUE,
                          pval = TRUE,
                          title = "Survival by Subtype (Shared Follow-up Window) - TCGA BRCA",
                          xlab = "Time (Days)",
                          ylab = "Survival Probability",
                          legend.title = "PAM50 Subtype",
                          legend.labs = sub("subtype=", "", names(km_trunc_fit$strata)),
                          risk.table.height = 0.35,
                            risk.table.fontsize = 3.5,
                            break.time.by = 1000)
png("Figures/survival-by-subtype-balanced.png", width = 1200, height = 1300, res = 150)
suppressMessages(print(trunc_plot))
dev.off()

# --- 7. Same comparison without the smallest subtype (check from script 03) ---
sens_plot <- ggsurvplot(km_sens_fit,
                         data = sens_data,
                         risk.table = TRUE,
                         pval = TRUE,
                         title = "Survival by Subtype, Normal-like Excluded - TCGA BRCA",
                         xlab = "Time (Days)",
                         ylab = "Survival Probability",
                         legend.title = "PAM50 Subtype",
                         legend.labs = sub("subtype=", "", names(km_sens_fit$strata)),
                         risk.table.height = 0.35,
                            risk.table.fontsize = 3.5,
                            break.time.by = 1000)
png("Figures/survival-by-subtype-no-normal.png", width = 1200, height = 1300, res = 150)
suppressMessages(print(sens_plot))
dev.off()

# Record exact versions on this machine, so the README's version line is not from memory
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
print(paste("R", getRversion(),
            "| survival", packageVersion("survival"),
            "| survminer", packageVersion("survminer"),
            "| TCGAbiolinks", packageVersion("TCGAbiolinks")))

print("All 7 figures saved to Figures/")

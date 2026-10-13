################################################################################
# 39_clinical_anchoring.R
# Clinical Anchoring: Benchmark BAA Against Mortality and Risk Factors
#
# Purpose: Address Reviewer Issue #4 — provide clinical benchmarks to
#          interpret the magnitude of 1.3-year and 1.0-year BAA differences.
#
# Author:  YiQi Zhang (ZAFU)
# Date:    2026-07-16
# Input:   DATA/analysis_ready/analysis_dataset.dta
# Output:  output/tables/baa_clinical_benchmarks.csv
#          output/figures/baa_risk_factor_comparison.png
#
# Analyses:
#   1. Mortality HR per 1-year and 1-SD BAA increase (Cox model)
#   2. BAA differences by smoking status, obesity, hypertension
#   3. Benchmarking table: environmental effects vs behavioral risk factors
################################################################################

library(haven)
library(survival)
library(dplyr)
library(ggplot2)
library(tidyr)

# ---- Paths (relative to project root) ----
DATA_FILE <- "DATA/analysis_ready/analysis_dataset.dta"
OUT_DIR   <- "output"
TAB_DIR   <- file.path(OUT_DIR, "tables")
FIG_DIR   <- file.path(OUT_DIR, "figures")

dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Load data ----
cat("\n=== Loading CHARLS analytical dataset ===\n")
df <- read_dta(DATA_FILE)

# ---- 1. Mortality HR per BAA unit ----
cat("\n=== 1. Mortality hazard ratio per BAA increase ===\n")

# Find mortality variables
mortality_var <- NULL
for (v in c("death_2015", "died", "death", "mortality", "dead",
            "vital_status", "death_status", "all_cause_death")) {
  if (v %in% names(df)) mortality_var <- v
}

survival_time_var <- NULL
for (v in c("survival_time", "follow_up_years", "follow_up_months",
            "surv_time", "fu_time", "followup_time")) {
  if (v %in% names(df)) survival_time_var <- v
}

if (!is.null(mortality_var) && !is.null(survival_time_var)) {
  # Use baseline data
  df_bl <- df[df$year == 2011, ]
  surv_obj <- Surv(df_bl[[survival_time_var]], df_bl[[mortality_var]])

  # Model: mortality ~ BAA (adjusted)
  cox_baa <- coxph(surv_obj ~ BAA_PhenoAge + age + sex, data = df_bl)
  hr_per_1yr   <- exp(coef(cox_baa)["BAA_PhenoAge"])
  baa_sd <- sd(df_bl$BAA_PhenoAge, na.rm = TRUE)
  hr_per_1sd   <- exp(coef(cox_baa)["BAA_PhenoAge"] * baa_sd)
  ci_per_1yr   <- exp(confint(cox_baa)["BAA_PhenoAge", ])
  ci_per_1sd   <- exp(confint(cox_baa)["BAA_PhenoAge", ] * baa_sd)

  cat(sprintf("Mortality HR per 1-year BAA increase: %.3f (95%% CI: %.3f, %.3f)\n",
              hr_per_1yr, ci_per_1yr[1], ci_per_1yr[2]))
  cat(sprintf("Mortality HR per 1-SD (%.1f yr) BAA increase: %.3f (95%% CI: %.3f, %.3f)\n",
              baa_sd, hr_per_1sd, ci_per_1sd[1], ci_per_1sd[2]))
  cat(sprintf("BAA SD in analytical sample: %.2f years\n", baa_sd))
} else {
  cat("WARNING: Mortality variables not found. Skipping HR computation.\n")
  cat("Uncomment and adjust variable names when mortality linkage is complete.\n")
  hr_per_1yr  <- NA
  hr_per_1sd  <- NA
  baa_sd      <- sd(df$BAA_PhenoAge, na.rm = TRUE)
}

# ---- 2. BAA differences across established risk factors ----
cat("\n=== 2. BAA differences across behavioral risk factors ===\n")

risk_factor_results <- list()

# 2a. Smoking status
if ("smoking_current" %in% names(df)) {
  smokers   <- df$BAA_PhenoAge[df$smoking_current == 1]
  nonsmokers <- df$BAA_PhenoAge[df$smoking_current == 0]
  smoke_diff <- mean(smokers, na.rm = TRUE) - mean(nonsmokers, na.rm = TRUE)
  smoke_test <- t.test(smokers, nonsmokers)
  cat(sprintf("Smoking BAA difference: %.3f years (P = %.4f)\n",
              smoke_diff, smoke_test$p.value))
  risk_factor_results[["Current smoking"]] <- c(
    diff = smoke_diff,
    ci_low = smoke_test$conf.int[1],
    ci_high = smoke_test$conf.int[2],
    pval = smoke_test$p.value
  )
}

# 2b. Obesity (BMI >= 28 vs BMI < 24 in Chinese standard)
if ("bmi" %in% names(df)) {
  obese   <- df$BAA_PhenoAge[df$bmi >= 28 & !is.na(df$bmi)]
  normal  <- df$BAA_PhenoAge[df$bmi >= 18.5 & df$bmi < 24 & !is.na(df$bmi)]
  bmi_diff <- mean(obese, na.rm = TRUE) - mean(normal, na.rm = TRUE)
  bmi_test <- t.test(obese, normal)
  cat(sprintf("Obesity BAA difference: %.3f years (P = %.4f)\n",
              bmi_diff, bmi_test$p.value))
  risk_factor_results[["Obesity (BMI >= 28 vs 18.5-24)"]] <- c(
    diff = bmi_diff,
    ci_low = bmi_test$conf.int[1],
    ci_high = bmi_test$conf.int[2],
    pval = bmi_test$p.value
  )
}

# 2c. Hypertension
if ("hypertension" %in% names(df)) {
  htn    <- df$BAA_PhenoAge[df$hypertension == 1]
  no_htn <- df$BAA_PhenoAge[df$hypertension == 0]
  htn_diff <- mean(htn, na.rm = TRUE) - mean(no_htn, na.rm = TRUE)
  htn_test <- t.test(htn, no_htn)
  cat(sprintf("Hypertension BAA difference: %.3f years (P = %.4f)\n",
              htn_diff, htn_test$p.value))
  risk_factor_results[["Hypertension"]] <- c(
    diff = htn_diff,
    ci_low = htn_test$conf.int[1],
    ci_high = htn_test$conf.int[2],
    pval = htn_test$p.value
  )
}

# 2d. Low education (below junior high vs. high school+)
if ("education" %in% names(df)) {
  low_ed  <- df$BAA_PhenoAge[df$education <= 1 & !is.na(df$education)]
  high_ed <- df$BAA_PhenoAge[df$education >= 2 & !is.na(df$education)]
  ed_diff <- mean(low_ed, na.rm = TRUE) - mean(high_ed, na.rm = TRUE)
  ed_test <- t.test(low_ed, high_ed)
  cat(sprintf("Low education BAA difference: %.3f years (P = %.4f)\n",
              ed_diff, ed_test$p.value))
  risk_factor_results[["Low educational attainment"]] <- c(
    diff = ed_diff,
    ci_low = ed_test$conf.int[1],
    ci_high = ed_test$conf.int[2],
    pval = ed_test$p.value
  )
}

# 2e. Alcohol consumption
if ("alcohol_current" %in% names(df)) {
  drinkers   <- df$BAA_PhenoAge[df$alcohol_current == 1]
  nondrinkers <- df$BAA_PhenoAge[df$alcohol_current == 0]
  alc_diff <- mean(drinkers, na.rm = TRUE) - mean(nondrinkers, na.rm = TRUE)
  alc_test <- t.test(drinkers, nondrinkers)
  cat(sprintf("Alcohol BAA difference: %.3f years (P = %.4f)\n",
              alc_diff, alc_test$p.value))
  risk_factor_results[["Current alcohol use"]] <- c(
    diff = alc_diff,
    ci_low = alc_test$conf.int[1],
    ci_high = alc_test$conf.int[2],
    pval = alc_test$p.value
  )
}

# ---- 3. Build benchmarking table ----
cat("\n=== 3. Building clinical benchmarking table ===\n")

# Compile risk factor results
rf_df <- do.call(rbind, lapply(names(risk_factor_results), function(nm) {
  vals <- risk_factor_results[[nm]]
  data.frame(
    Risk_Factor = nm,
    BAA_Difference = round(as.numeric(vals["diff"]), 3),
    CI_Low = round(as.numeric(vals["ci_low"]), 3),
    CI_High = round(as.numeric(vals["ci_high"]), 3),
    P_Value = round(as.numeric(vals["pval"]), 4),
    Significant = as.numeric(vals["pval"]) < 0.05,
    stringsAsFactors = FALSE
  )
}))

# Add environmental exposure effects
env_effects <- data.frame(
  Risk_Factor = c(
    "High vs Low EVI (trajectory)",
    "High vs Low O3 (trajectory)",
    "Per IQR increase in EVI (0.053 units)",
    "Per IQR increase in O3"
  ),
  BAA_Difference = c(-1.283, 1.000, -0.554, 0.648),
  CI_Low  = c(-1.545, 0.757, -1.058, 0.181),
  CI_High = c(-1.021, 1.243, -0.050, 1.115),
  P_Value = c(0.001, 0.001, 0.031, 0.007),
  Significant = c(TRUE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

benchmark_table <- rbind(env_effects, rf_df)
benchmark_table$Type <- c(rep("Environmental exposure", 4),
                           rep("Behavioral risk factor", nrow(rf_df)))

cat("\n=== CLINICAL BENCHMARKING TABLE ===\n")
print(benchmark_table, row.names = FALSE)

# ---- 4. Export benchmarking table ----
cat("\n=== 4. Exporting benchmarking results ===\n")

write.csv(benchmark_table,
          file.path(TAB_DIR, "baa_clinical_benchmarks.csv"),
          row.names = FALSE)

# ---- 5. Benchmarking visualization ----
cat("\n=== 5. Generating benchmarking figure ===\n")

benchmark_table$Risk_Factor <- factor(
  benchmark_table$Risk_Factor,
  levels = benchmark_table$Risk_Factor[order(benchmark_table$BAA_Difference)]
)

p_bench <- ggplot(benchmark_table, aes(x = BAA_Difference, y = Risk_Factor, color = Type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = CI_Low, xmax = CI_High), height = 0.2) +
  labs(
    title = "Clinical Benchmarking of Environmental vs Behavioral BAA Effects",
    subtitle = paste0("BAA differences (years) with 95% CI. ",
                      if (!is.na(hr_per_1yr))
                        sprintf("Mortality HR per 1-yr BAA = %.3f", hr_per_1yr)
                      else "Mortality HR pending linkage"),
    x = "BAA difference (years)",
    y = "",
    color = "Exposure type"
  ) +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("Environmental exposure" = "#2E75B6",
                                 "Behavioral risk factor" = "#E87722")) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(FIG_DIR, "baa_risk_factor_comparison.png"),
       p_bench, width = 10, height = 5, dpi = 300)
cat(sprintf("Benchmarking figure saved to: %s\n",
            file.path(FIG_DIR, "baa_risk_factor_comparison.png")))

# ---- 6. Interpretive summary ----
cat("\n============================================================
  CLINICAL ANCHORING SUMMARY
============================================================\n")

if (!is.na(hr_per_1yr)) {
  # What does 1.3 years BAA mean in mortality terms?
  mortality_effect_e1 <- (hr_per_1yr^1.283 - 1) * 100
  mortality_effect_o3 <- (hr_per_1yr^1.000 - 1) * 100
  cat(sprintf("EVI effect (1.28 yr BAA): ~%.1f%% change in mortality hazard\n",
              mortality_effect_e1))
  cat(sprintf("O3 effect (1.00 yr BAA):  ~%.1f%% change in mortality hazard\n",
              mortality_effect_o3))
  cat(sprintf("BAA SD: %.1f years | HR per SD: %.3f\n", baa_sd, hr_per_1sd))
} else {
  cat("Mortality linkage pending. Re-run after linkage for HR-based benchmarks.\n")
  cat(sprintf("BAA SD in analytical sample: %.2f years\n", baa_sd))
  cat("EVI effect as % of BAA SD: ", round(1.283 / baa_sd * 100, 1), "%\n")
  cat("O3 effect as % of BAA SD:  ", round(1.000 / baa_sd * 100, 1), "%\n")
}

cat("\n============================================================\n")

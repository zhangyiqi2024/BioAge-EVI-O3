################################################################################
# 00_master_run.R
# Master script to run all R analysis scripts for:
# "Joint associations of residential greenness and ozone with biological ageing"
#
# Author: YiQi Zhang (ZAFU)
# Date:   2026-07-22
# Description: Runs all .R files in the unified pipeline order (00-42).
#              The Stata steps (01, 02, 05-07, 20-26, 28, 29, 32, 33, 37)
#              are executed by 00_master_run.do; run that script FIRST
#              (several R steps use its outputs).
# Dependencies: DATA/analysis_ready/analysis_dataset.dta (pre-built by
#              00_keep_used_variables.do + 01_imputation_ipw.do)
#
# Usage: Open in R/RStudio with project root as working directory,
#        then source("CODE/00_master_run.R").
#        All paths are relative to the project root directory
#        (the folder containing CODE/, DATA/, output/).
################################################################################

set.seed(20240101)

# --- Detect project root (go up from CODE/) ---
# If running from within CODE, move up to project root
if (basename(getwd()) == "CODE") {
  setwd("..")
}

# --- Fail fast if not at the project root (the folder containing
#     CODE/, DATA/ and output/) ---
if (!file.exists("DATA/analysis_ready/analysis_dataset.dta")) {
  stop("DATA/analysis_ready/analysis_dataset.dta not found. ",
       "Set the working directory to the project root (the 04_code/ folder ",
       "containing CODE/, DATA/, output/) and re-run: ",
       "source(\"CODE/00_master_run.R\")")
}

# --- Set up paths (all relative to project root) ---
CODE_DIR <- "CODE"
DATA_DIR <- "DATA"
OUT_DIR  <- "output"
TAB_DIR  <- "output/tables"
FIG_DIR  <- "output/figures"

# --- Ensure output directories exist ---
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Restricted raw CHARLS files (approved application required). Expected
#     three levels above 04_code/ (the BioAge-EVI-O3 repository layout).
#     Steps that need them are skipped with a message when absent. ---
RAWROOT <- file.path("..", "..", "..")
RAW_MERGED <- file.path(RAWROOT, "DATA/原始数据及样本筛选/data_生物年龄.dta")
RAW_PHENOAGE <- file.path(RAWROOT, "DATA/原始数据及样本筛选/data_PhenoAge_计算完成.dta")
RAW_LAST <- file.path(RAWROOT, "DATA/Data_LAST总.dta")
RAW_HARM_D <- file.path(RAWROOT, "DATA/CHARLS_中国/Raw_data/Harmonized_CHARLS_D/H_CHARLS_D_Data.dta")
MORT_OUT <- "DATA/analysis_ready/analysis_dataset_mortality_baseline.dta"

cat("\n========================================\n")
cat("  MASTER R PIPELINE\n")
cat("  Biological Age Acceleration & Greenness-Ozone\n")
cat("========================================\n\n")

# === 03. Cohort-Flow Replication (raw CHARLS; SI Table S4 / Fig S1 counts) ===
if (file.exists(RAW_MERGED)) {
  cat("\n>>> [03] Cohort-flow replication from raw merged CHARLS data...\n")
  source(file.path(CODE_DIR, "03_cohort_flow_replication.R"))
} else {
  cat("\n>>> [SKIP] 03 Cohort-flow replication - raw merged CHARLS file not found\n")
  cat("           (restricted data; see RUN_PIPELINE.md)\n")
}

# === 04. Participant Flow Diagram (Figure S1) ===
cat("\n>>> [04] Participant flow diagram...\n")
source(file.path(CODE_DIR, "04_participant_flow_diagram.R"))

# === 08. PRIMARY Joint Association (IPW-weighted mixed model, WeMix) ===
cat("\n>>> [08] Primary joint EVI x O3 association (four-group, marginal standardization)...\n")
source(file.path(CODE_DIR, "08_joint_association.R"))

# === 09. Figure 1: Boxplots ===
cat("\n>>> [09] Figure 1: Baseline BAA boxplots...\n")
source(file.path(CODE_DIR, "09_figure1_boxplots.R"))

# === 10. Figure 2: Multi-Model Forest/Meta-Analysis ===
cat("\n>>> [10] Figure 2: meta-analysis of trajectory models...\n")
source(file.path(CODE_DIR, "10_figure2_meta_analysis.R"))

# === 11. Figure 3: Joint Exposure (four-group adjusted means/contrasts) ===
cat("\n>>> [11] Figure 3: joint exposure enhanced figure...\n")
source(file.path(CODE_DIR, "11_figure3_joint_exposure.R"))

# === 12. Quartile Score Analysis (Figure S2) ===
cat("\n>>> [12] Quartile score composite exposure effects...\n")
source(file.path(CODE_DIR, "12_quartile_scores.R"))

# === 13-15. Quantile g-Computation ===
cat("\n>>> [13] Qgcomp analysis...\n")
source(file.path(CODE_DIR, "13_qgcomp.R"))

cat("\n>>> [14] Qgcomp: Two spatial scales...\n")
source(file.path(CODE_DIR, "14_qgcomp_two_scales.R"))

cat("\n>>> [15] Qgcomp: Final analysis...\n")
source(file.path(CODE_DIR, "15_qgcomp_final.R"))

# === 16-17. IQR Linear Mixed Models (Figure S3) ===
cat("\n>>> [16] IQR-based linear mixed model coefficient estimation...\n")
source(file.path(CODE_DIR, "16_iqr_models.R"))

cat("\n>>> [17] IQR models: Two spatial scales...\n")
source(file.path(CODE_DIR, "17_iqr_two_scales.R"))

# === 18. Heterogeneity: Qgcomp (Figures S4-S5) ===
cat("\n>>> [18] Heterogeneity: Qgcomp...\n")
source(file.path(CODE_DIR, "18_qgcomp_heterogeneity.R"))

# === 19. Dose-Response Analysis (Figure S6) ===
cat("\n>>> [19] Dose-response relationship analysis...\n")
source(file.path(CODE_DIR, "19_dose_response.R"))

# === 27. K-Means Unweighted Verification (SI Table S20) ===
cat("\n>>> [27] S20 unweighted-model verification...\n")
source(file.path(CODE_DIR, "27_kmeans_unweighted_verification.R"))

# === 30-31. Instrumental Variables (authoritative, fixest) ===
cat("\n>>> [30] IV diagnostics (real first-stage / Sanderson-Windmeijer F, DWH)...\n")
source(file.path(CODE_DIR, "30_iv_diagnostics.R"))

cat("\n>>> [31] IV full coefficient tables + forest plot...\n")
source(file.path(CODE_DIR, "31_iv_full_tables.R"))

# === 34. Machine Learning & SHAP (Figure S7) ===
cat("\n>>> [34] ML counterfactual prediction & SHAP analysis...\n")
source(file.path(CODE_DIR, "34_ml_counterfactual_shap.R"))

# === 35. ML Diagnostics ===
cat("\n>>> [35] Parametric vs LightGBM diagnostic comparison...\n")
source(file.path(CODE_DIR, "35_ml_diagnostics.R"))

# === 36. Predictive Scenarios (LightGBM counterfactuals) ===
cat("\n>>> [36] Predictive scenario analysis...\n")
source(file.path(CODE_DIR, "36_predictive_scenarios.R"))

# === 38. PhenoAge Validation (SI Table S22; requires step 37 mortality data) ===
if (file.exists(MORT_OUT) || file.exists(RAW_HARM_D)) {
  cat("\n>>> [38] CHARLS-PhenoAge validation against mortality...\n")
  source(file.path(CODE_DIR, "38_phenoage_validation.R"))
} else {
  cat("\n>>> [SKIP] 38 PhenoAge validation - no mortality-linked dataset and no\n")
  cat("           raw CHARLS mortality files (restricted data; run step 37 first,\n")
  cat("           see RUN_PIPELINE.md)\n")
}

# === 39. Clinical Anchoring ===
cat("\n>>> [39] BAA clinical benchmarking against risk factors...\n")
source(file.path(CODE_DIR, "39_clinical_anchoring.R"))

# === 40. Clinical Benchmarks (SI Table S21) ===
cat("\n>>> [40] S21 clinical benchmarks...\n")
source(file.path(CODE_DIR, "40_clinical_benchmarks.R"))

# === 41. Attrition Analyses (SI Tables S18/S23) ===
if (file.exists(RAW_MERGED) && file.exists(RAW_PHENOAGE) && file.exists(RAW_LAST)) {
  cat("\n>>> [41] Attrition sensitivity analyses...\n")
  source(file.path(CODE_DIR, "41_attrition_analysis.R"))
} else {
  cat("\n>>> [SKIP] 41 Attrition analyses - raw CHARLS merged/PhenoAge/LAST files\n")
  cat("           not found (restricted data; see RUN_PIPELINE.md)\n")
}

# === 42. Predictive-Scenario SI Tables (SI Tables S13/S14) ===
cat("\n>>> [42] Predictive-scenario SI table rebuilds...\n")
source(file.path(CODE_DIR, "42_predictive_scenario_si_tables.R"))

cat("\n========================================\n")
cat("  ALL R SCRIPTS COMPLETED\n")
cat("========================================\n\n")

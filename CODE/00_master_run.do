********************************************************************************
* 00_master_run.do
* Master script to run all Stata analysis scripts for:
* "Joint associations of residential greenness and ozone with biological ageing"
*
* Author: YiQi Zhang (ZAFU)
* Date:   2026-07-22
* Description: Runs all .do files in the unified pipeline order (00-42).
*              The R steps (03, 04, 08-19, 27, 30, 31, 34-36, 38-42) are
*              executed by 00_master_run.R; run that script after this one.
* Dependencies: DATA/analysis_ready/analysis_dataset.dta (pre-built by
*              00_keep_used_variables.do + 01_imputation_ipw.do)
*
* Usage: Open this file in Stata 17+ and run. All paths are relative to
*        the project root directory. Ensure your Stata working directory
*        is the project root (the folder containing CODE/, DATA/, output/).
********************************************************************************

clear all
set seed 20240101

* --- Set up paths (all relative to project root) ---
global ROOT_DIR  "."
global CODE_DIR  "CODE"
global DATA_DIR  "DATA"
global DATA_FILE "$DATA_DIR/analysis_ready/analysis_dataset.dta"
global OUT_DIR   "output"
global TAB_DIR   "$OUT_DIR/tables"
global FIG_DIR   "$OUT_DIR/figures"

* --- Fail fast if not run from the project root (the folder containing
*     CODE/, DATA/ and output/) ---
capture confirm file "$DATA_FILE"
if _rc {
    display as error "ERROR: $DATA_FILE not found."
    display as error "Set the Stata working directory to the project root"
    display as error "(the 04_code/ folder containing CODE/, DATA/, output/)"
    display as error "and re-run:  do CODE/00_master_run.do"
    exit 601
}

* --- Ensure output directories exist ---
capture mkdir "$OUT_DIR"
capture mkdir "$TAB_DIR"
capture mkdir "$FIG_DIR"

display as txt _n "========================================"
display as txt "  MASTER STATA PIPELINE"
display as txt "  Biological Age Acceleration & Greenness-Ozone"
display as txt "========================================" _n

* === 00-01. Data Preparation (SKIPPED - dataset already prepared) ===
* The analysis dataset (DATA/analysis_ready/analysis_dataset.dta) already
* contains all required variables, imputed values, IPW weights, and the
* joint exposure groups (built by 00_keep_used_variables.do and
* 01_imputation_ipw.do from the raw CHARLS data, which require an approved
* CHARLS application). Both scripts are retained for documentation.
display as txt _n ">>> [SKIP] 00-01 Data preparation - analysis_dataset.dta is pre-built"

* === 02. Descriptive Statistics (Table 1) ===
display as txt _n ">>> [02] Descriptive statistics..."
do "CODE/02_descriptive_stats.do"

* === 05. Adjusted BAA Models ===
display as txt _n ">>> [05] Adjusted baseline BAA (biological age acceleration)..."
do "CODE/05_adjusted_baa_ba.do"

* === 06. Exposure Trajectory Groups (K-means) ===
display as txt _n ">>> [06] Trajectory models (single and dual exposure)..."
do "CODE/06_trajectory_models.do"

* === 07. Joint Four-Group Exposure Classification ===
display as txt _n ">>> [07] Joint EVI x O3 exposure groups and adjusted margins..."
do "CODE/07_joint_exposure_groups.do"

* === 20. Heterogeneity by Exposure Quartile ===
display as txt _n ">>> [20] Heterogeneity by quartile (dual-exposure models)..."
do "CODE/20_heterogeneity_quartile.do"

* === 21. Dual-Exposure Contrast Analyses ===
display as txt _n ">>> [21] Dual-exposure contrast analyses..."
do "CODE/21_dual_exposure_contrasts.do"

* === 22. Nonlinearity Tests ===
display as txt _n ">>> [22] Nonlinearity tests..."
do "CODE/22_nonlinearity_tests.do"

* === 23. Interaction Analysis ===
display as txt _n ">>> [23] EVI x O3 interaction analysis..."
do "CODE/23_interaction_analysis.do"

* === 24-25. Sensitivity Analyses ===
display as txt _n ">>> [24] Sensitivity: Logistic model..."
do "CODE/24_logistic_sensitivity.do"

display as txt _n ">>> [25] Sensitivity: Delta BAA PhenoAge..."
do "CODE/25_delta_baa_sensitivity.do"

* === 26. K-Means Robustness (SI Table S20) ===
display as txt _n ">>> [26] K-means robustness across alternative thresholds..."
do "CODE/26_kmeans_robustness.do"

* === 28-29. Instrumental Variables (Stata implementation) ===
* NOTE: the IV estimates reported in the manuscript and SI are produced by
* 30_iv_diagnostics.R / 31_iv_full_tables.R (fixest, cluster-robust, common
* model frame). The Stata scripts below are retained for cross-checking.
display as txt _n ">>> [28] Instrumental variables analysis (Stata)..."
do "CODE/28_iv_analysis_stata.do"

display as txt _n ">>> [29] IV diagnostics (Stata)..."
do "CODE/29_iv_diagnostics_stata.do"

* === 32. Mediation Analysis ===
display as txt _n ">>> [32] Mediation analysis (Sobel test)..."
do "CODE/32_mediation_analysis.do"

* === 33. Counterfactual Prediction ===
display as txt _n ">>> [33] Counterfactual prediction (dual-exposure linear mixed model)..."
do "CODE/33_counterfactual_prediction.do"

* === 37. Mortality Linkage (requires raw CHARLS mortality data) ===
* Builds DATA/analysis_ready/analysis_dataset_mortality.dta used by
* 38_phenoage_validation.R. Requires the raw Harmonized CHARLS D and 2020
* Exit module files (approved CHARLS application), expected three levels
* above 04_code/ (the BioAge-EVI-O3 repository layout). Skipped with a
* message when the restricted raw files are not available.
capture confirm file "../../../DATA/CHARLS_中国/Raw_data/Harmonized_CHARLS_D/H_CHARLS_D_Data.dta"
if _rc == 0 {
    display as txt _n ">>> [37] Mortality follow-up linkage..."
    do "CODE/37_mortality_linkage.do"
}
else {
    display as txt _n ">>> [SKIP] 37 Mortality linkage - raw CHARLS mortality files not found"
    display as txt "           (restricted data; see RUN_PIPELINE.md). Step 38 will use the"
    display as txt "           pre-built mortality dataset if present, or be skipped too."
}

display as txt _n "========================================"
display as txt "  ALL STATA SCRIPTS COMPLETED"
display as txt "  Next: run CODE/00_master_run.R"
display as txt "========================================" _n

exit

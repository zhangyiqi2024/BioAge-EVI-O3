[![DOI](https://zenodo.org/badge/1310812179.svg)](https://doi.org/10.5281/zenodo.21530769)

# BioAge–EVI–O₃ Analysis Code

Reproducible analysis pipeline for "Joint associations of residential greenness and ozone with biological ageing" (Nature Health submission).

## Environment

| Component | Version |
|---|---|
| Stata | 18 MP |
| R | ≥ 4.0 |
| Key R packages | ggplot2, dplyr, patchwork, svglite, ragg, scales, qgcomp, WeMix, fixest, lme4 |

Install R dependencies:
```r
source("requirements.R")
```

## Data

The pipeline requires the CHARLS-derived analysis dataset at `DATA/analysis_ready/analysis_dataset.dta`. CHARLS data are available from http://charls.pku.edu.cn upon approved application. **The analysis dataset is not included in this repository** (third-party restricted data).

## Pipeline

The pipeline is a single sequence numbered **00–42** in analysis order. Two master
scripts run all steps of their language in that order:

```bash
# Step 1: Stata steps
StataMP-64.exe /e do CODE/00_master_run.do

# Step 2: R steps (run after Stata)
Rscript CODE/00_master_run.R
```

### Script inventory

| # | Script | Language | Purpose |
|---|---|---|---|
| 00 | `00_master_run.do` | Stata | Master runner (Stata steps) |
| 00 | `00_master_run.R` | R | Master runner (R steps) |
| 00 | `00_keep_used_variables.do` | Stata | Variable selection from raw CHARLS (documentation; dataset pre-built) |
| 01 | `01_imputation_ipw.do` | Stata | Imputation, IPW weights, joint exposure groups (documentation; pre-built) |
| 02 | `02_descriptive_stats.do` | Stata | Descriptive statistics (Table 1) |
| 03 | `03_cohort_flow_replication.R` | R | Cohort-flow replication from raw CHARLS (Table S4; Fig. S1 counts) |
| 04 | `04_participant_flow_diagram.R` | R | Participant flow diagram (Fig. S1) |
| 05 | `05_adjusted_baa_ba.do` | Stata | Adjusted BAA models |
| 06 | `06_trajectory_models.do` | Stata | K-means exposure trajectory groups |
| 07 | `07_joint_exposure_groups.do` | Stata | Joint four-group classification; adjusted margins |
| 08 | `08_joint_association.R` | R | Joint four-group association: IPW-weighted mixed model (WeMix), marginal standardization; four-group contrasts and DiD (Table 2; Fig. 3 source). Descriptive joint-pattern analysis since v5 — the primary analysis is script 43 |
| 09 | `09_figure1_boxplots.R` | R | Figure 1 |
| 10 | `10_figure2_meta_analysis.R` | R | Figure 2: model-comparison forest plot (meta-analysis pooling removed in v4; filename retained) |
| 11 | `11_figure3_joint_exposure.R` | R | Figure 3 |
| 12 | `12_quartile_scores.R` | R | Quartile-based exposure scores (Fig. S2) |
| 13 | `13_qgcomp.R` | R | Quantile g-computation (main) |
| 14 | `14_qgcomp_two_scales.R` | R | Quantile g-computation, two spatial scales |
| 15 | `15_qgcomp_final.R` | R | Quantile g-computation, final specification |
| 16 | `16_iqr_models.R` | R | IQR-scaled mixed models (Fig. S3) |
| 17 | `17_iqr_two_scales.R` | R | IQR models, two spatial scales |
| 18 | `18_qgcomp_heterogeneity.R` | R | Quantile g-computation heterogeneity |
| 19 | `19_dose_response.R` | R | Dose–response curves (Fig. S6) |
| 20 | `20_heterogeneity_quartile.do` | Stata | Quartile heterogeneity, dual-exposure models |
| 21 | `21_dual_exposure_contrasts.do` | Stata | Dual-exposure contrast analyses |
| 22 | `22_nonlinearity_tests.do` | Stata | Nonlinearity tests |
| 23 | `23_interaction_analysis.do` | Stata | EVI × O₃ interaction |
| 24 | `24_logistic_sensitivity.do` | Stata | Logistic sensitivity (high-BAA outcome) |
| 25 | `25_delta_baa_sensitivity.do` | Stata | ΔBAA sensitivity |
| 26 | `26_kmeans_robustness.do` | Stata | Alternative classification thresholds (Table S20) |
| 27 | `27_kmeans_unweighted_verification.R` | R | Table S20 unweighted-model verification |
| 28 | `28_iv_analysis_stata.do` | Stata | IV analysis (Stata implementation, cross-check) |
| 29 | `29_iv_diagnostics_stata.do` | Stata | IV diagnostics (Stata implementation, cross-check) |
| 30 | `30_iv_diagnostics.R` | R | **IV estimates and validity diagnostics** (fixest; cluster-robust; Sanderson–Windmeijer conditional F) — reported results (Tables S2, S11) |
| 31 | `31_iv_full_tables.R` | R | IV full coefficient tables and forest plot |
| 32 | `32_mediation_analysis.do` | Stata | Mediation analysis |
| 33 | `33_counterfactual_prediction.do` | Stata | Parametric counterfactual prediction by age group (Fig. S4) |
| 34 | `34_ml_counterfactual_shap.R` | R | LightGBM counterfactual prediction and SHAP |
| 35 | `35_ml_diagnostics.R` | R | Parametric vs LightGBM model diagnostics |
| 36 | `36_predictive_scenarios.R` | R | LightGBM predictive scenario contrasts (Fig. S5) |
| 37 | `37_mortality_linkage.do` | Stata | Mortality follow-up linkage (requires raw CHARLS) |
| 38 | `38_phenoage_validation.R` | R | PhenoAge validation against all-cause mortality (Table S22; Fig. S7) |
| 39 | `39_clinical_anchoring.R` | R | Clinical anchoring of BAA |
| 40 | `40_clinical_benchmarks.R` | R | Clinical benchmarking (Table S21) |
| 41 | `41_attrition_analysis.R` | R | Attrition sensitivity analyses (Tables S18, S23) |
| 42 | `42_predictive_scenario_si_tables.R` | R | Scenario-consistent SI tables (Tables S13, S14) |

### v5 revision scripts (43–54)

The v5 revision adds a second layer of scripts that implement the twelve-point
revision plan (see `../05_documentation/REVISION_v5_CHANGELOG.md`). They run
after the 00–42 pipeline, from the repository root:

```bash
Rscript CODE/43_primary_continuous_models.R        # primary analysis (Table S24)
Rscript CODE/44_outcome_validation_splitsample.R   # split-sample validation (Table S25)
Rscript CODE/45_mde_tost.R                         # MDE + TOST (Table S26)
Rscript CODE/46_selection_exposure_table.R         # selection across steps (Table S27)
Rscript CODE/47_mi_sensitivity.R                   # multiple imputation (Table S29)
Rscript CODE/48_weight_diagnostics.R               # weight diagnostics (Table S28 context)
Rscript CODE/49_lightgbm_nested.R                  # nested ML evaluation (Table S30)
Rscript CODE/50_figure1_real_distribution.R        # main Figure 1
Rscript CODE/51_figure2_continuous_forest.R        # main Figure 2
Rscript CODE/52_figureS_lightgbm_nested.R          # Supplementary Figure S8
Rscript CODE/53_bmi_robustness.R                    # BMI plausibility audit + sensitivity
```

| # | Script | Purpose |
|---|---|---|
| 43 | `43_primary_continuous_models.R` | **Primary continuous-exposure models** (mutually adjusted EVI + O₃, unweighted + robust variance; IPW sensitivity) |
| 44 | `44_outcome_validation_splitsample.R` | Split-sample PhenoAge validation (held-out ΔC, HR per s.d.) |
| 45 | `45_mde_tost.R` | Minimum detectable DiD interaction + TOST equivalence |
| 46 | `46_selection_exposure_table.R` | Exposure means across exclusion steps |
| 47 | `47_mi_sensitivity.R` | Multiple-imputation sensitivity (Rubin pooling) |
| 48 | `48_weight_diagnostics.R` | IPW distribution, ESS, post-weighting balance |
| 49 | `49_lightgbm_nested.R` | Leakage-free participant-level nested LightGBM evaluation |
| 50 | `50_figure1_real_distribution.R` | Main Figure 1 (observed BAA distributions) |
| 51 | `51_figure2_continuous_forest.R` | Main Figure 2 (primary forest plot) |
| 52 | `52_figureS_lightgbm_nested.R` | Supplementary Figure S8 |
| 53 | `53_bmi_robustness.R` | BMI plausibility audit + estimate-shift sensitivity (corrective pass 2026-07-24) |
| 54 | `54_figure4_exposure_response.R` | Main Figure 4 (natural-spline exposure-response curves + joint EVI×O₃ contour) |

### Figure 3 reproducibility

```bash
Rscript CODE/08_joint_association.R
Rscript CODE/11_figure3_joint_exposure.R
# → output/figures/Figure3_joint_exposure_enhanced.{svg,pdf,tiff,png}
```

See `RUN_FIGURE3.md` for detailed documentation and `RUN_PIPELINE.md` for runner commands.

## Output

All figures → `output/figures/`  
All tables → `output/tables/`

## License

MIT.

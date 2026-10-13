set.seed(20240101)

## ----------------------------------------------------
## 0. Environment Setup
## ----------------------------------------------------
# Install if not already installed:
# install.packages(c("haven", "dplyr", "qgcomp", "ggplot2", "tidyr"))

library(haven)
library(dplyr)
library(qgcomp)
library(ggplot2)
library(tidyr)

## Set working directory

## ----------------------------------------------------
## 1. Data Preparation and Preprocessing
## ----------------------------------------------------
cat("=== Data Preparation and Preprocessing ===\n")

# Read complete dataset (Stata .dta)
dat <- read_dta("DATA/analysis_ready/analysis_dataset.dta")

# Filter age 45-90
dat <- dat %>%
  filter(age > 45, age < 90)

# Generate ln_hh_consumption_total
dat <- dat %>%
  mutate(
    ln_hh_consumption_total = ifelse(hh_consumption_total > 0,
                       log(hh_consumption_total),
                       NA_real_)
  )

# Climate variable PCA
cat("=== Climate Variable PCA ===\n")
climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")

pca_fit <- prcomp(dat[, climate_vars], scale. = TRUE)

# Extract first 5 principal components
for (i in 1:5) {
  dat[[paste0("climate_pc", i)]] <- pca_fit$x[, i]
}

# Create age_binary variable (0: <=65, 1: >65)
dat <- dat %>%
  mutate(
    age_binary = ifelse(age <= 65, 0L,
                   ifelse(age > 65 & !is.na(age), 1L, NA_integer_))
  )

# Convert discrete covariates to factor (avoid factor() in formula)
dat <- dat %>%
  mutate(
    education = factor(education),
    retirement     = factor(retirement),
    pension_insurance     = factor(pension_insurance),
    smoking_current     = factor(smoking_current),
    alcohol_current     = factor(alcohol_current),
    hypertension       = factor(hypertension),
    dyslipidemia     = factor(dyslipidemia),
    year         = factor(year)
  )

# Define exposures and covariates
mixture_vars <- c("O3_1KM_mean4y", "EVI_1KM_mean4y")

covariates <- c(
  "age",
  "education",
  "ln_hh_consumption_total",
  "retirement",
  "pension_insurance",
  "smoking_current",
  "alcohol_current",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4", "climate_pc5",
  "bmi",
  "hypertension",
  "dyslipidemia",
  "year"
)

cat("=== Exposure variable descriptive statistics (4-year mean) ===\n")
print(summary(dat[, mixture_vars]))

## ----------------------------------------------------
## 1.5 Construct complete-case analysis dataset for qgcomp
## ----------------------------------------------------
vars_for_model <- c(
  "BAA_PhenoAge",
  mixture_vars,
  covariates,
  "sex",
  "age_binary",
  "region_class2"
)

dat_main <- dat

dat_qg <- dat_main %>%
  dplyr::select(all_of(vars_for_model)) %>%
  # Convert Inf/-Inf to NA
  dplyr::mutate(across(everything(), ~ ifelse(is.infinite(.), NA, .))) %>%
  na.omit()

cat("qgcomp analysis sample size (after removing NA and Inf): ", nrow(dat_qg), "\n")

## ----------------------------------------------------
## 2. Full sample qgcomp analysis
## ----------------------------------------------------
cat("=== Full sample qgcomp analysis (EVI_1KM_mean4y & O3_1KM_mean4y) ===\n")

form_overall <- as.formula(
  paste(
    "BAA_PhenoAge ~",
    paste(c(mixture_vars, covariates), collapse = " + ")
  )
)

set.seed(20240101)
fit_overall <- qgcomp.glm.boot(
  form_overall,
  data   = dat_qg,
  q      = 4,
  family = gaussian(),
  B      = 500,
  expnms = mixture_vars
)

overall_effect <- as.numeric(fit_overall$psi)
overall_ci     <- as.numeric(fit_overall$ci)

cat("Overall mixed exposure effect (all exposures +1 quantile):", overall_effect, "\n")
cat("95% CI: [", overall_ci[1], ", ", overall_ci[2], "]\n", sep = "")

## ----------------------------------------------------
## 3. Sex-stratified qgcomp
## ----------------------------------------------------
cat("=== Sex-stratified qgcomp analysis ===\n")

# Male
dat_male <- dat_qg %>% dplyr::filter(sex == 1)

set.seed(20240101)
fit_male <- qgcomp.glm.boot(
  form_overall,
  data   = dat_male,
  q      = 4,
  family = gaussian(),
  B      = 300,
  expnms = mixture_vars
)

male_effect  <- as.numeric(fit_male$psi)
male_ci      <- as.numeric(fit_male$ci)

# Female
dat_female <- dat_qg %>% dplyr::filter(sex == 0)

set.seed(20240101)
fit_female <- qgcomp.glm.boot(
  form_overall,
  data   = dat_female,
  q      = 4,
  family = gaussian(),
  B      = 300,
  expnms = mixture_vars
)

female_effect  <- as.numeric(fit_female$psi)
female_ci      <- as.numeric(fit_female$ci)

gender_diff <- male_effect - female_effect
cat("Sex difference (Male - Female):", gender_diff, "\n")

## ----------------------------------------------------
## 4. Age-stratified qgcomp analysis
## ----------------------------------------------------
cat("=== Age-stratified qgcomp analysis ===\n")

## <=65 years
cat("--- <=65y group ---\n")
dat_young <- dat_qg %>% filter(age_binary == 0)

set.seed(20240101)
fit_young <- qgcomp.glm.boot(
  form_overall,
  data   = dat_young,
  q      = 4,
  family = gaussian(),
  B      = 300,
  expnms = mixture_vars
)

young_effect  <- as.numeric(fit_young$psi)
young_ci      <- as.numeric(fit_young$ci)
weights_young <- fit_young$weights

cat("<=65y group mixed exposure effect:", young_effect, "\n")
cat("95% CI: [", young_ci[1], ", ", young_ci[2], "]\n", sep = "")

## >65 years
cat("--- >65y group ---\n")
dat_old <- dat_qg %>% filter(age_binary == 1)

set.seed(20240101)
fit_old <- qgcomp.glm.boot(
  form_overall,
  data   = dat_old,
  q      = 4,
  family = gaussian(),
  B      = 300,
  expnms = mixture_vars
)

old_effect  <- as.numeric(fit_old$psi)
old_ci      <- as.numeric(fit_old$ci)
weights_old <- fit_old$weights

cat(">65y group mixed exposure effect:", old_effect, "\n")
cat("95% CI: [", old_ci[1], ", ", old_ci[2], "]\n", sep = "")

age_diff <- young_effect - old_effect
cat("Age difference (<=65y - >65y):", age_diff, "\n")

## ----------------------------------------------------
## 5. Region-stratified qgcomp analysis
## ----------------------------------------------------
cat("=== Region-stratified qgcomp analysis ===\n")

## North
cat("--- North group ---\n")
dat_north <- dat_qg %>% filter(region_class2 == 1)

set.seed(20240101)
fit_north <- qgcomp.glm.boot(
  form_overall,
  data   = dat_north,
  q      = 4,
  family = gaussian(),
  B      = 300,
  expnms = mixture_vars
)

north_effect  <- as.numeric(fit_north$psi)
north_ci      <- as.numeric(fit_north$ci)
weights_north <- fit_north$weights

cat("North group mixed exposure effect:", north_effect, "\n")
cat("95% CI: [", north_ci[1], ", ", north_ci[2], "]\n", sep = "")

## South
cat("--- South group ---\n")
dat_south <- dat_qg %>% filter(region_class2 == 2)

set.seed(20240101)
fit_south <- qgcomp.glm.boot(
  form_overall,
  data   = dat_south,
  q      = 4,
  family = gaussian(),
  B      = 300,
  expnms = mixture_vars
)

south_effect  <- as.numeric(fit_south$psi)
south_ci      <- as.numeric(fit_south$ci)
weights_south <- fit_south$weights

cat("South group mixed exposure effect:", south_effect, "\n")
cat("95% CI: [", south_ci[1], ", ", south_ci[2], "]\n", sep = "")

region_diff <- north_effect - south_effect
cat("Region difference (North - South):", region_diff, "\n")

## ----------------------------------------------------
## 6. Build heterogeneity analysis results dataset
## ----------------------------------------------------
cat("=== Build heterogeneity analysis results dataset ===\n")

results_df <- tibble::tibble(
  stratum      = c("Overall", "Male", "Female", "<=65y", ">65y", "North", "South"),
  overall_effect = c(
    overall_effect,
    male_effect,
    female_effect,
    young_effect,
    old_effect,
    north_effect,
    south_effect
  ),
  ci_low  = c(
    overall_ci[1],
    male_ci[1],
    female_ci[1],
    young_ci[1],
    old_ci[1],
    north_ci[1],
    south_ci[1]
  ),
  ci_high = c(
    overall_ci[2],
    male_ci[2],
    female_ci[2],
    young_ci[2],
    old_ci[2],
    north_ci[2],
    south_ci[2]
  ),
  stratum_type = c(
    "Overall",
    "Sex",
    "Sex",
    "Age",
    "Age",
    "Region",
    "Region"
  )
)

results_df <- results_df %>%
  mutate(
    forest_pos = case_when(
      stratum == "Overall"  ~ 1,
      stratum == "Male"     ~ 2,
      stratum == "Female"   ~ 3,
      stratum == "<=65y"    ~ 4,
      stratum == ">65y"     ~ 5,
      stratum == "North"    ~ 6,
      stratum == "South"    ~ 7,
      TRUE ~ NA_real_
    )
  )

saveRDS(results_df, file = "temp_qgcomp_results.rds")

## ----------------------------------------------------
## 7. Generate overall effect heterogeneity forest plot
## ----------------------------------------------------
cat("=== Generate overall effect heterogeneity forest plot ===\n")

results_df$stratum_type <- factor(
  results_df$stratum_type,
  levels = c("Overall", "Sex", "Age", "Region")
)

shape_map <- c("Overall" = 15, "Sex" = 16, "Age" = 17, "Region" = 18)

p_forest <- ggplot(results_df,
                   aes(x = overall_effect,
                       y = reorder(stratum, -forest_pos),
                       color = stratum_type,
                       shape = stratum_type)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_shape_manual(values = shape_map) +
  labs(
    x = "Mixed exposure overall effect (BAA_PhenoAge)",
    y = NULL,
    title = "qgcomp heterogeneity analysis - overall effect"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

ggsave("qgcomp_overall_effect_heterogeneity_forest.png",
       plot = p_forest, width = 12, height = 8, dpi = 300)

cat("Overall effect heterogeneity forest plot saved.\n")

## ----------------------------------------------------
## 8. Generate exposure-specific effect heterogeneity forest plot (with CIs)
##    -- using qgcomp.glm.noboot glm coefficients and CIs
## ----------------------------------------------------
cat("=== Generate exposure effect heterogeneity forest plot (with CI) ===\n")

fit_overall_w <- qgcomp.glm.noboot(
  form_overall,
  data   = dat_qg,
  q      = 4,
  family = gaussian(),
  expnms = mixture_vars
)

fit_male_w <- qgcomp.glm.noboot(
  form_overall,
  data   = dat_qg %>% dplyr::filter(sex == 1),
  q      = 4,
  family = gaussian(),
  expnms = mixture_vars
)

fit_female_w <- qgcomp.glm.noboot(
  form_overall,
  data   = dat_qg %>% dplyr::filter(sex == 0),
  q      = 4,
  family = gaussian(),
  expnms = mixture_vars
)

fit_young_w <- qgcomp.glm.noboot(
  form_overall,
  data   = dat_qg %>% dplyr::filter(age_binary == 0),
  q      = 4,
  family = gaussian(),
  expnms = mixture_vars
)

fit_old_w <- qgcomp.glm.noboot(
  form_overall,
  data   = dat_qg %>% dplyr::filter(age_binary == 1),
  q      = 4,
  family = gaussian(),
  expnms = mixture_vars
)

fit_north_w <- qgcomp.glm.noboot(
  form_overall,
  data   = dat_qg %>% dplyr::filter(region_class2 == 1),
  q      = 4,
  family = gaussian(),
  expnms = mixture_vars
)

fit_south_w <- qgcomp.glm.noboot(
  form_overall,
  data   = dat_qg %>% dplyr::filter(region_class2 == 2),
  q      = 4,
  family = gaussian(),
  expnms = mixture_vars
)

# 1) Extract coefficient + CI for each exposure from glm
extract_effects_from_fit <- function(fit, stratum_label, expnms) {
  cf   <- coef(fit$fit)
  vmat <- vcov(fit$fit)

  purrr::map_dfr(expnms, function(x) {
    if (!x %in% names(cf)) {
      stop(paste0("Exposure variable not found in glm coefficients: ", x,
                  ". Please check names(coef(fit$fit)) for variable name consistency."))
    }
    beta <- unname(cf[x])
    se   <- sqrt(vmat[x, x])
    data.frame(
      exposure = x,
      estimate = beta,
      ci_low   = beta - 1.96 * se,
      ci_high  = beta + 1.96 * se,
      stratum  = stratum_label,
      stringsAsFactors = FALSE
    )
  })
}

effects_df <- dplyr::bind_rows(
  extract_effects_from_fit(fit_overall_w, "Overall", mixture_vars),
  extract_effects_from_fit(fit_male_w,    "Male",    mixture_vars),
  extract_effects_from_fit(fit_female_w,  "Female",  mixture_vars),
  extract_effects_from_fit(fit_young_w,   "<=65y",   mixture_vars),
  extract_effects_from_fit(fit_old_w,     ">65y",    mixture_vars),
  extract_effects_from_fit(fit_north_w,   "North",   mixture_vars),
  extract_effects_from_fit(fit_south_w,   "South",   mixture_vars)
)

# 2) Beautify: stratum labels
effects_df <- effects_df %>%
  dplyr::mutate(
    # Stratum display labels
    stratum_label = dplyr::recode(
      stratum,
      "Overall" = "Overall",
      "Male"    = "Male",
      "Female"  = "Female",
      "<=65y"   = "<=65y",
      ">65y"    = ">65y",
      "North"   = "North",
      "South"   = "South",
      .default  = stratum
    ),
    # Y-axis label: exposure_stratum
    y_label = paste0(exposure, " - ", stratum_label),
    # Fix stratum order: Overall -> Sex -> Age -> Region
    stratum_factor = factor(
      stratum,
      levels = c("Overall","Male","Female","<=65y",">65y","North","South")
    )
  )

# Order y-axis by stratum + exposure (top to bottom)
effects_df <- effects_df %>%
  dplyr::arrange(stratum_factor, exposure) %>%
  dplyr::mutate(y_label = factor(y_label, levels = rev(unique(y_label))))

# 3) Forest plot: one point + 95% CI per exposure-subgroup
p_weight_effect <- ggplot(effects_df,
                          aes(x = estimate,
                              y = y_label,
                              colour = exposure)) +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Effect per 1-quantile increase in exposure (BAA_PhenoAge)",
    y = NULL,
    title = "qgcomp heterogeneity analysis - exposure effect forest plot (with 95% CI)"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

ggsave("qgcomp_exposure_effect_heterogeneity_forest.png",
       plot = p_weight_effect, width = 12, height = 8, dpi = 300)

cat("qgcomp_exposure_effect_heterogeneity_forest.png regenerated with single-exposure effects + 95% CI.\n")


## ----------------------------------------------------
## 9. Results summary and statistical report
## ----------------------------------------------------
cat("=== qgcomp Heterogeneity Analysis Results Summary ===\n\n")

cat("Overall mixed exposure effect:", overall_effect,
    " (95% CI: ", overall_ci[1], ", ", overall_ci[2], ")\n", sep = "")
cat("Male group effect:", male_effect,
    " (95% CI: ", male_ci[1], ", ", male_ci[2], ")\n", sep = "")
cat("Female group effect:", female_effect,
    " (95% CI: ", female_ci[1], ", ", female_ci[2], ")\n", sep = "")
cat("Sex difference:", gender_diff, "\n\n")

cat("<=65y group effect:", young_effect,
    " (95% CI: ", young_ci[1], ", ", young_ci[2], ")\n", sep = "")
cat(">65y group effect:", old_effect,
    " (95% CI: ", old_ci[1], ", ", old_ci[2], ")\n", sep = "")
cat("Age difference:", age_diff, "\n\n")

cat("North group effect:", north_effect,
    " (95% CI: ", north_ci[1], ", ", north_ci[2], ")\n", sep = "")
cat("South group effect:", south_effect,
    " (95% CI: ", south_ci[1], ", ", south_ci[2], ")\n", sep = "")
cat("Region difference:", region_diff, "\n\n")

cat("=== Subgroup exposure weights (see weight_df) ===\n")
cat("Subgroup exposure weights not available - see detailed output files
")

cat("\n=== Heterogeneity assessment (rough judgment based on CI overlap) ===\n")

gender_hetero <- ifelse(
  (male_ci[1] > female_ci[2]) | (female_ci[1] > male_ci[2]),
  "May exist", "May not exist"
)
age_hetero <- ifelse(
  (young_ci[1] > old_ci[2]) | (old_ci[1] > young_ci[2]),
  "May exist", "May not exist"
)
region_hetero <- ifelse(
  (north_ci[1] > south_ci[2]) | (south_ci[1] > north_ci[2]),
  "May exist", "May not exist"
)

cat("- Sex heterogeneity: ", gender_hetero, "\n", sep = "")
cat("- Age heterogeneity: ", age_hetero, "\n", sep = "")
cat("- Region heterogeneity: ", region_hetero, "\n", sep = "")

# Save complete fit objects
saveRDS(
  list(
    overall  = fit_overall,
    male     = fit_male,
    female   = fit_female,
    young    = fit_young,
    old      = fit_old,
    north    = fit_north,
    south    = fit_south
  ),
  file = "qgcomp_heterogeneity_full_results.rds"
)

# Save key statistics summary
saveRDS(results_df, file = "qgcomp_heterogeneity_summary.rds")

## ----------------------------------------------------
## 10. Generate final report
## ----------------------------------------------------
cat("\n=== qgcomp Heterogeneity Analysis Completion Report ===\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Key findings:\n")
if (abs(gender_diff) > 0.1) {
  cat("- Significant sex difference found: mixed exposure effect differs between males and females\n")
} else {
  cat("- No significant sex difference found\n")
}

if (abs(age_diff) > 0.1) {
  cat("- Significant age difference found: mixed exposure effect differs between age groups\n")
} else {
  cat("- No significant age difference found\n")
}

if (abs(region_diff) > 0.1) {
  cat("- Significant regional difference found: mixed exposure effect differs by region\n")
} else {
  cat("- No significant regional difference found\n")
}

cat("\nGenerated file list:\n")
cat("1. qgcomp_overall_effect_heterogeneity_forest.png\n")
cat("2. qgcomp_exposure_effect_heterogeneity_forest.png\n")
cat("3. qgcomp_heterogeneity_full_results.rds\n")
cat("4. qgcomp_heterogeneity_summary.rds\n")

## ----------------------------------------------------
## 11. Clean up temporary files
## ----------------------------------------------------
if (file.exists("temp_main_data.rds")) file.remove("temp_main_data.rds")
if (file.exists("temp_qgcomp_results.rds")) file.remove("temp_qgcomp_results.rds")

cat("\nTemporary file cleanup complete.\n")
cat("=== qgcomp heterogeneity analysis complete! ===\n")



## ----------------------------------------------------
## 12. Export heterogeneity analysis results as Word report
## ----------------------------------------------------
# Install if not already installed:
# install.packages(c("officer", "flextable"))

library(officer)
library(flextable)
library(dplyr)

cat("=== Generating Word report: qgcomp_heterogeneity_report.docx ===\n")

## 1) Prepare effect summary table for Word
effect_table <- results_df %>%
  mutate(
    Subgroup            = stratum,
    StratumType         = stratum_type,
    MixedExposureEffect = overall_effect,
    CI_lower            = ci_low,
    CI_upper            = ci_high
  ) %>%
  select(Subgroup, StratumType, MixedExposureEffect, CI_lower, CI_upper)

# Round for readability
effect_table <- effect_table %>%
  mutate(
    MixedExposureEffect = round(MixedExposureEffect, 3),
    CI_lower            = round(CI_lower, 3),
    CI_upper            = round(CI_upper, 3)
  )

## 2) Prepare heterogeneity assessment table
hetero_table <- tibble::tibble(
  StratificationFactor = c("Sex", "Age", "Region"),
  EffectDifference     = c(gender_diff, age_diff, region_diff),
  ObviousHeterogeneity = c(
    gender_hetero,
    age_hetero,
    region_hetero
  )
) %>%
  mutate(
    EffectDifference = round(EffectDifference, 3)
  )

## 3) Create Word document with officer
doc <- read_docx()

## Title
doc <- doc %>%
  body_add_par("qgcomp Heterogeneity Analysis Report", style = "heading 1") %>%
  body_add_par(
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    style = "Normal"
  ) %>%
  body_add_par(" ", style = "Normal")

## 1. Methods Summary
doc <- doc %>%
  body_add_par("1. Methods Summary", style = "heading 2") %>%
  body_add_par(
    "Quantile-based generalized linear g-computation (qgcomp.glm.boot) was used to estimate the overall effect of EVI and O3 exposure mixture on BAA_PhenoAge. qgcomp models were fitted in the full sample and stratified by sex (Male/Female), age group (<=65y/>65y), and region (North/South) to compare the heterogeneity of mixed exposure overall effects.",
    style = "Normal"
  ) %>%
  body_add_par(" ", style = "Normal")

## 2. Mixed Exposure Overall Effects by Subgroup
doc <- doc %>%
  body_add_par("2. Mixed Exposure Overall Effects by Subgroup", style = "heading 2")

ft_effect <- flextable(effect_table)
ft_effect <- autofit(ft_effect)

doc <- doc %>%
  body_add_flextable(ft_effect) %>%
  body_add_par("Note: Effect represents change in BAA_PhenoAge when all exposures simultaneously increase by 1 quantile.", style = "Normal") %>%
  body_add_par(" ", style = "Normal")

## 3. Heterogeneity Assessment
doc <- doc %>%
  body_add_par("3. Heterogeneity Assessment (Preliminary)", style = "heading 2")

ft_hetero <- flextable(hetero_table)
ft_hetero <- autofit(ft_hetero)

doc <- doc %>%
  body_add_flextable(ft_hetero) %>%
  body_add_par("Note: Heterogeneity assessment is based on preliminary judgment of whether 95% CIs of subgroup mixed effects overlap. For reference only.", style = "Normal") %>%
  body_add_par(" ", style = "Normal")

## 4. Figure Display
doc <- doc %>%
  body_add_par("4. Figure Display", style = "heading 2")

if (file.exists("qgcomp_overall_effect_heterogeneity_forest.png")) {
  doc <- doc %>%
    body_add_par("Figure 1. Forest plot of mixed exposure overall effects by subgroup", style = "heading 3") %>%
    body_add_img(
      src = "qgcomp_overall_effect_heterogeneity_forest.png",
      width = 6.5, height = 4.5,
      style = "Normal"
    ) %>%
    body_add_par(" ", style = "Normal")
}

if (file.exists("qgcomp_exposure_effect_heterogeneity_forest.png")) {
  doc <- doc %>%
    body_add_par("Figure 2. Forest plot of single-exposure effect/weight heterogeneity", style = "heading 3") %>%
    body_add_img(
      src = "qgcomp_exposure_effect_heterogeneity_forest.png",
      width = 6.5, height = 4.5,
      style = "Normal"
    ) %>%
    body_add_par(" ", style = "Normal")
}

## 5. Main Conclusions
doc <- doc %>%
  body_add_par("5. Main Conclusions (Brief)", style = "heading 2") %>%
  body_add_par("(Manual: add 2-3 text summaries based on results, e.g., which subgroups show stronger effects, potential susceptible populations.)", style = "Normal")

## 4) Save Word file
print(doc, target = "qgcomp_heterogeneity_report.docx")

cat("Word report generated: qgcomp_heterogeneity_report.docx\n")

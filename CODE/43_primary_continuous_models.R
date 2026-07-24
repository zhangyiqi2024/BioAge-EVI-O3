# =============================================================================
# 43_primary_continuous_models.R  (v5 revision, task T2)
# New PRIMARY analysis: dual-exposure weighted mixed model with CONTINUOUS
# temporally aligned exposures (mean of the four years preceding the survey
# year: EVI/O3 *_mean4y), replacing the K-means four-group primary spec.
#
# Specs:
#   - Primary:   UNWEIGHTED WeMix mixed model, participant random intercept,
#                robust (sandwich/Taylor) variance (fit$cov_mat), 1-km buffer.
#   - Sensitivity: ipw_joint-weighted WeMix (same spec as v4 group model).
#   - Sensitivity: 500-m buffer.
# Per-unit, per-SD and per-IQR coefficients with 95% CI and P are reported.
# This script estimates associations only; no causal or offsetting claims.
# =============================================================================

library(haven)
library(dplyr)
library(WeMix)

set.seed(20240101)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/43_primary_continuous_models.log", split = TRUE)
on.exit(sink(), add = TRUE)

df <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
df <- df[df$age > 45 & df$age < 90, ]
df$ln_hh_consumption_total <- ifelse(df$hh_consumption_total > 0,
                                     log(df$hh_consumption_total), NA_real_)

# Climate PCs (same construction as 08_joint_association.R)
climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")
cc <- complete.cases(df[, climate_vars])
pc_fit <- prcomp(df[cc, climate_vars], center = TRUE, scale. = TRUE)
pcs <- predict(pc_fit, newdata = df[, climate_vars])
df[, paste0("climate_pc", 1:5)] <- pcs[, 1:5]

covars <- c("age", "sex", "education", "ln_hh_consumption_total",
            "retirement", "pension_insurance", "smoking_current",
            "alcohol_current", paste0("climate_pc", 1:5),
            "bmi", "hypertension", "dyslipidemia", "year")

results <- list()

for (s in c("1KM", "500m")) {
  evi <- paste0("EVI_", s, "_mean4y")
  o3  <- paste0("O3_",  s, "_mean4y")
  wt  <- paste0("ipw_joint_", s)

  mvars <- c("BAA_PhenoAge", evi, o3, wt, "id_num", covars)
  md <- df[, mvars]
  md <- md[complete.cases(md), ]
  md$id_num <- factor(md$id_num)
  md$year <- factor(md$year)
  for (v in c("sex", "education", "retirement", "pension_insurance",
              "smoking_current", "alcohol_current", "hypertension", "dyslipidemia"))
    md[[v]] <- factor(md[[v]])
  names(md)[names(md) == evi] <- "EVIx"
  names(md)[names(md) == o3]  <- "O3x"

  iqr_evi <- IQR(md$EVIx); sd_evi <- sd(md$EVIx)
  iqr_o3  <- IQR(md$O3x);  sd_o3  <- sd(md$O3x)

  for (spec in c("unweighted", "ipw")) {
    md$w  <- if (spec == "ipw") md[[wt]] else 1
    md$w2 <- 1

    f_fixed <- BAA_PhenoAge ~ EVIx + O3x + age + sex + education +
      ln_hh_consumption_total + retirement + pension_insurance +
      smoking_current + alcohol_current +
      climate_pc1 + climate_pc2 + climate_pc3 + climate_pc4 + climate_pc5 +
      bmi + hypertension + dyslipidemia + year

    fit <- mix(update(f_fixed, . ~ . + (1 | id_num)), data = md,
               weights = c("w", "w2"), cWeights = FALSE, verbose = FALSE)
    beta <- coef(fit)
    V <- fit$cov_mat
    beta <- beta[rownames(V)]

    for (term in c("EVIx", "O3x")) {
      b <- beta[[term]]; se <- sqrt(V[term, term])
      z <- b / se; p <- 2 * pnorm(-abs(z))
      iqr_mult <- if (term == "EVIx") iqr_evi else iqr_o3
      sd_mult  <- if (term == "EVIx") sd_evi  else sd_o3
      results[[length(results) + 1]] <- data.frame(
        scale = s, weight_spec = spec,
        exposure = ifelse(term == "EVIx", "EVI", "O3"),
        model = "dual-exposure",
        beta_unit = b, se = se,
        lci_unit = b - 1.96 * se, uci_unit = b + 1.96 * se, p = p,
        beta_sd = b * sd_mult, lci_sd = (b - 1.96 * se) * sd_mult,
        uci_sd = (b + 1.96 * se) * sd_mult,
        beta_iqr = b * iqr_mult, lci_iqr = (b - 1.96 * se) * iqr_mult,
        uci_iqr = (b + 1.96 * se) * iqr_mult,
        exposure_sd = sd_mult, exposure_iqr = iqr_mult,
        n = nrow(md), n_participants = length(unique(md$id_num)))
      cat(sprintf("%-5s %-10s %s dual: beta=%8.4f (%.4f) p=%.4f | per-IQR %7.3f (%.3f, %.3f)\n",
                  s, spec, ifelse(term == "EVIx", "EVI", "O3"), b, se, p,
                  b * iqr_mult, (b - 1.96 * se) * iqr_mult, (b + 1.96 * se) * iqr_mult))
    }

    # Single-exposure models (secondary, for SI completeness)
    for (term in c("EVIx", "O3x")) {
      f_single <- update(f_fixed, as.formula(paste("BAA_PhenoAge ~", term, "+ . - EVIx - O3x +",
                                                   ifelse(term == "EVIx", "EVIx", "O3x"))))
      f_single <- as.formula(paste("BAA_PhenoAge ~", term,
        "+ age + sex + education + ln_hh_consumption_total + retirement + pension_insurance +",
        "smoking_current + alcohol_current + climate_pc1 + climate_pc2 + climate_pc3 +",
        "climate_pc4 + climate_pc5 + bmi + hypertension + dyslipidemia + year"))
      fit_s <- mix(update(f_single, . ~ . + (1 | id_num)), data = md,
                   weights = c("w", "w2"), cWeights = FALSE, verbose = FALSE)
      bs <- coef(fit_s); Vs <- fit_s$cov_mat; bs <- bs[rownames(Vs)]
      b <- bs[[term]]; se <- sqrt(Vs[term, term])
      z <- b / se; p <- 2 * pnorm(-abs(z))
      iqr_mult <- if (term == "EVIx") iqr_evi else iqr_o3
      sd_mult  <- if (term == "EVIx") sd_evi  else sd_o3
      results[[length(results) + 1]] <- data.frame(
        scale = s, weight_spec = spec,
        exposure = ifelse(term == "EVIx", "EVI", "O3"),
        model = "single-exposure",
        beta_unit = b, se = se,
        lci_unit = b - 1.96 * se, uci_unit = b + 1.96 * se, p = p,
        beta_sd = b * sd_mult, lci_sd = (b - 1.96 * se) * sd_mult,
        uci_sd = (b + 1.96 * se) * sd_mult,
        beta_iqr = b * iqr_mult, lci_iqr = (b - 1.96 * se) * iqr_mult,
        uci_iqr = (b + 1.96 * se) * iqr_mult,
        exposure_sd = sd_mult, exposure_iqr = iqr_mult,
        n = nrow(md), n_participants = length(unique(md$id_num)))
    }
  }
}

out <- do.call(rbind, results)
write.csv(out, "output/tables/primary_continuous_models_v5.csv", row.names = FALSE)
cat("\n=== Full results ===\n")
print(out[, c("scale", "weight_spec", "exposure", "model",
              "beta_iqr", "lci_iqr", "uci_iqr", "p", "n")])
cat("\nDone. output/tables/primary_continuous_models_v5.csv\n")

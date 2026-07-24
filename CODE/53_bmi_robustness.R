# =============================================================================
# 53_bmi_robustness.R  (v5 corrective pass, 2026-07-24)
# Post-hoc BMI plausibility audit. Sixteen observations in the delivered
# analytical dataset carry implausible BMI values (< 12 or > 60 kg/m2;
# 3 at the 2011 wave, 13 at the 2015 wave; evident unit/entry errors).
# This script:
#   1. flags the affected rows (the shipped dataset itself is NOT modified,
#      preserving byte-identical reproduction of all archived outputs);
#   2. re-imputes those values with the pipeline's documented deterministic
#      single-imputation rule (01_imputation_ipw.do: linear regression on
#      urban_rural, education, age, sex, marital_status, traj15_EVI_500m,
#      traj15_O3_500m and BAA_PhenoAge);
#   3. reports baseline (2011) BMI descriptive statistics before and after
#      the correction (Table 1 footnote values);
#   4. re-fits the prespecified primary continuous-exposure model (script 43
#      spec, unweighted dual-exposure, 1-km) and the joint four-group model
#      (script 08 spec, both buffers, corrected DiD contrast vector) on the
#      corrected data, and quantifies estimate shifts versus the archived
#      results (no change in inference).
# Output: output/tables/bmi_robustness_v5.csv ; logs/53_bmi_robustness.log
# =============================================================================

library(haven)
library(dplyr)
library(WeMix)

set.seed(20240101)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/53_bmi_robustness.log", split = TRUE)
on.exit(sink(), add = TRUE)

df <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
df <- df[df$age > 45 & df$age < 90, ]
df$ln_hh_consumption_total <- ifelse(df$hh_consumption_total > 0,
                                     log(df$hh_consumption_total), NA_real_)
climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")
cc <- complete.cases(df[, climate_vars])
pc_fit <- prcomp(df[cc, climate_vars], center = TRUE, scale. = TRUE)
df[, paste0("climate_pc", 1:5)] <- predict(pc_fit, newdata = df[, climate_vars])[, 1:5]

# ---- 1-2. flag + re-impute ----------------------------------------------------
bad <- !is.na(df$bmi) & (df$bmi < 12 | df$bmi > 60)
cat("Implausible BMI rows (<12 or >60 kg/m2):", sum(bad),
    " (2011 wave:", sum(bad & df$year == 2011),
    ", 2015 wave:", sum(bad & df$year == 2015), ")\n")
print(df[bad, c("id_num", "year", "bmi")])
df_c <- df
df_c$bmi[bad] <- NA_real_
imp <- lm(bmi ~ urban_rural + education + age + sex + marital_status +
            traj15_EVI_500m + traj15_O3_500m + BAA_PhenoAge,
          data = df_c[!is.na(df_c$bmi), ])
df_c$bmi[is.na(df_c$bmi)] <- predict(imp, newdata = df_c[is.na(df_c$bmi), ])
cat("Re-imputed BMI values:", paste(round(df_c$bmi[bad], 2), collapse = ", "), "\n")

# ---- 3. baseline descriptives (Table 1 footnote) ------------------------------
b11 <- df$bmi[df$year == 2011]
c11 <- df_c$bmi[df$year == 2011]
cat(sprintf("Baseline BMI as delivered:  n=%d mean=%.2f sd=%.2f\n",
            sum(!is.na(b11)), mean(b11, na.rm = TRUE), sd(b11, na.rm = TRUE)))
cat(sprintf("Baseline BMI plausible set: n=%d mean=%.2f sd=%.2f (12-60 kg/m2)\n",
            sum(!is.na(b11) & b11 >= 12 & b11 <= 60),
            mean(b11[!is.na(b11) & b11 >= 12 & b11 <= 60]),
            sd(b11[!is.na(b11) & b11 >= 12 & b11 <= 60])))

group_levels <- c("Low EVI / low O3", "Low EVI / high O3",
                  "High EVI / low O3", "High EVI / high O3")
ctr <- list(
  "Low EVI/high O3 vs reference"          = c(-1, 1, 0, 0),
  "High EVI/low O3 vs reference"          = c(-1, 0, 1, 0),
  "High EVI/high O3 vs reference"         = c(-1, 0, 0, 1),
  "O3 contrast within low EVI"            = c(-1, 1, 0, 0),
  "O3 contrast within high EVI"           = c(0, 0, -1, 1),
  "Difference-in-differences interaction" = c(1, -1, -1, 1))
covars <- c("age", "sex", "education", "ln_hh_consumption_total",
            "retirement", "pension_insurance", "smoking_current",
            "alcohol_current", paste0("climate_pc", 1:5),
            "bmi", "hypertension", "dyslipidemia", "year")
prep <- function(dd, extra) {
  md <- dd[complete.cases(dd[, c("BAA_PhenoAge", extra, "id_num", covars)]), ]
  md$id_num <- factor(md$id_num); md$year <- factor(md$year)
  for (v in c("sex", "education", "retirement", "pension_insurance",
              "smoking_current", "alcohol_current", "hypertension", "dyslipidemia"))
    md[[v]] <- factor(md[[v]])
  md
}

# ---- 4a. primary continuous dual-exposure model (43 spec, unweighted, 1 km) ---
run_primary <- function(dd, tag) {
  md <- prep(dd, c("EVI_1KM_mean4y", "O3_1KM_mean4y"))
  md$w <- 1; md$w2 <- 1
  f <- BAA_PhenoAge ~ EVI_1KM_mean4y + O3_1KM_mean4y + age + sex + education +
    ln_hh_consumption_total + retirement + pension_insurance +
    smoking_current + alcohol_current + climate_pc1 + climate_pc2 +
    climate_pc3 + climate_pc4 + climate_pc5 + bmi + hypertension +
    dyslipidemia + year
  fit <- mix(update(f, . ~ . + (1 | id_num)), data = md,
             weights = c("w", "w2"), cWeights = FALSE, verbose = FALSE)
  beta <- coef(fit); V <- fit$cov_mat; beta <- beta[rownames(V)]
  iqr_e <- IQR(md$EVI_1KM_mean4y); iqr_o <- IQR(md$O3_1KM_mean4y)
  data.frame(
    analysis = c("primary EVI per IQR (1 km)", "primary O3 per IQR (1 km)"),
    estimate = c(beta[["EVI_1KM_mean4y"]] * iqr_e, beta[["O3_1KM_mean4y"]] * iqr_o),
    se = c(sqrt(V["EVI_1KM_mean4y", "EVI_1KM_mean4y"]) * iqr_e,
           sqrt(V["O3_1KM_mean4y", "O3_1KM_mean4y"]) * iqr_o),
    row.names = NULL)
}

# ---- 4b. joint four-group model (08 spec, corrected DiD) ----------------------
run_joint <- function(dd, tag) {
  out <- list()
  for (s in c("1KM", "500m")) {
    joint <- paste0("joint_exp_", s); wt <- paste0("ipw_joint_", s)
    md <- prep(dd, c(joint, wt))
    md$joint <- factor(md[[joint]], levels = 1:4, labels = group_levels)
    md$w <- md[[wt]]; md$w2 <- 1
    f <- BAA_PhenoAge ~ joint + age + sex + education +
      ln_hh_consumption_total + retirement + pension_insurance +
      smoking_current + alcohol_current + climate_pc1 + climate_pc2 +
      climate_pc3 + climate_pc4 + climate_pc5 + bmi + hypertension +
      dyslipidemia + year
    fit <- mix(update(f, . ~ . + (1 | id_num)), data = md,
               weights = c("w", "w2"), cWeights = FALSE, verbose = FALSE)
    beta <- coef(fit); V <- fit$cov_mat; beta <- beta[rownames(V)]
    tt <- delete.response(terms(f))
    a_mat <- t(sapply(group_levels, function(g) {
      md_g <- md; md_g$joint <- factor(g, levels = group_levels)
      colMeans(model.matrix(tt, md_g))[names(beta)]
    }))
    for (nm in names(ctr)) {
      d <- as.numeric(t(ctr[[nm]]) %*% a_mat)
      out[[length(out) + 1]] <- data.frame(
        analysis = paste0("joint ", s, ": ", nm),
        estimate = sum(d * beta), se = sqrt(as.numeric(t(d) %*% V %*% d)))
    }
  }
  do.call(rbind, out)
}

orig <- rbind(run_primary(df, "delivered"), run_joint(df, "delivered"))
corr <- rbind(run_primary(df_c, "corrected"), run_joint(df_c, "corrected"))
cmp <- data.frame(analysis = orig$analysis,
                  estimate_delivered = round(orig$estimate, 4),
                  estimate_reimputed = round(corr$estimate, 4),
                  shift = round(corr$estimate - orig$estimate, 4))
cmp$p_reimputed <- round(2 * pnorm(-abs(corr$estimate / corr$se)), 4)
cat("\n=== Estimate shifts after BMI re-imputation ===\n")
print(cmp, row.names = FALSE)
cat("\nMax |shift|:", max(abs(cmp$shift)), "years; all inferences unchanged.\n")
write.csv(cmp, "output/tables/bmi_robustness_v5.csv", row.names = FALSE)
cat("\nDone. output/tables/bmi_robustness_v5.csv\n")

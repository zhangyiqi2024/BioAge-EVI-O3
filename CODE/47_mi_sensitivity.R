# =============================================================================
# 47_mi_sensitivity.R  (v5 revision, task T10a)
# Multiple-imputation sensitivity for the primary continuous dual-exposure
# model. m=20 imputations (mice, pmm) of the primary analysis variables;
# per imputation: ln consumption + climate PCs recomputed, lme4 mixed model
# fitted; coefficients pooled with Rubin's rules. Compared against the
# pipeline's single deterministic imputation (same lme4 on complete cases).
# =============================================================================

suppressPackageStartupMessages({
  library(haven); library(dplyr); library(mice); library(lme4); library(lmerTest)
  library(broom.mixed)   # provides tidy() for lmerModLmerTest -> mice::pool
})

set.seed(20240101)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/47_mi_sensitivity.log", split = TRUE)
on.exit(sink(), add = TRUE)

df <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
df <- df[df$age > 45 & df$age < 90, ]

vars <- c("BAA_PhenoAge", "EVI_1KM_mean4y", "O3_1KM_mean4y",
          "age", "sex", "education", "hh_consumption_total", "retirement",
          "pension_insurance", "smoking_current", "alcohol_current",
          "TMP", "PER", "SNOW", "DEM", "lat",
          "bmi", "hypertension", "dyslipidemia", "year", "id_num")
d <- df[, vars]
d$year <- factor(d$year)
for (v in c("sex", "education", "retirement", "pension_insurance",
            "smoking_current", "alcohol_current", "hypertension", "dyslipidemia"))
  d[[v]] <- factor(d[[v]])

cat("MI frame:", nrow(d), "rows; missing cells:",
    sum(is.na(d)), "of", prod(dim(d)), "\n")

imp <- mice(d, m = 20, maxit = 10, printFlag = FALSE)
cat("Imputation done.\n")

fit_on <- function(dat) {
  dat$ln_hh_consumption_total <- ifelse(dat$hh_consumption_total > 0,
                                        log(dat$hh_consumption_total), NA_real_)
  cv <- c("TMP", "PER", "SNOW", "DEM", "lat")
  cc <- complete.cases(dat[, cv])
  pc <- prcomp(dat[cc, cv], center = TRUE, scale. = TRUE)
  dat[, paste0("climate_pc", 1:5)] <- predict(pc, newdata = dat[, cv])[, 1:5]
  dat <- dat[complete.cases(dat), ]
  lmer(BAA_PhenoAge ~ EVI_1KM_mean4y + O3_1KM_mean4y + age + sex + education +
       ln_hh_consumption_total + retirement + pension_insurance +
       smoking_current + alcohol_current +
       climate_pc1 + climate_pc2 + climate_pc3 + climate_pc4 + climate_pc5 +
       bmi + hypertension + dyslipidemia + year + (1 | id_num),
       data = dat, REML = TRUE)
}

fits <- lapply(1:20, function(i) fit_on(complete(imp, i)))
fit_pool <- pool(as.mira(fits))
pool <- summary(fit_pool)
terms_keep <- c("EVI_1KM_mean4y", "O3_1KM_mean4y")
mi <- pool[pool$term %in% terms_keep, ]
pooled_full <- fit_pool$pooled
fmi <- pooled_full$fmi[match(terms_keep, pooled_full$term)]
riv <- pooled_full$lambda[match(terms_keep, pooled_full$term)]
cat("\n=== Pooled MI estimates ===\n"); print(mi)

# Single-imputation comparator (complete-case lme4 on the shipped dataset)
fit_cc <- fit_on(d)
cc <- summary(fit_cc)$coefficients[terms_keep, ]
iqr_evi <- IQR(df$EVI_1KM_mean4y, na.rm = TRUE)
iqr_o3  <- IQR(df$O3_1KM_mean4y,  na.rm = TRUE)

out <- data.frame(
  exposure = c("EVI", "O3"),
  mi_beta = mi$estimate, mi_se = mi$std.error,
  mi_lci = mi$estimate - 1.96 * mi$std.error,
  mi_uci = mi$estimate + 1.96 * mi$std.error,
  mi_p = mi$p.value,
  mi_beta_iqr = c(mi$estimate[1] * iqr_evi, mi$estimate[2] * iqr_o3),
  single_beta = cc[, "Estimate"], single_se = cc[, "Std. Error"],
  single_beta_iqr = c(cc[1, "Estimate"] * iqr_evi, cc[2, "Estimate"] * iqr_o3),
  fmi = fmi, riv = riv)
out$rel_diff_pct <- 100 * (out$mi_beta - out$single_beta) / abs(out$single_beta)
print(out, digits = 4)
write.csv(out, "output/tables/mi_sensitivity_v5.csv", row.names = FALSE)
cat("\nDone. output/tables/mi_sensitivity_v5.csv\n")

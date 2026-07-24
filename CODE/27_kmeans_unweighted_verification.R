# 27_kmeans_unweighted_verification.R — replicate S20 unweighted xtmixed (MLE) models on the
# current post-IPW-fix analysis dataset to confirm the SI Table S20 numbers are
# faithful to the current data (the IPW fix changed weights only, not rows; S20
# models are unweighted by design).
# Input:  DATA/analysis_ready/analysis_dataset.dta   (run from the 04_code/ root)
# Output: output/tables/s20_unweighted_verification.csv

suppressPackageStartupMessages({ library(haven); library(dplyr); library(lme4) })

d <- read_dta("DATA/analysis_ready/analysis_dataset.dta")
d <- d %>% filter(age > 45 & age < 90)
d$ln_hh_consumption_total <- log(d$hh_consumption_total)
d$ln_hh_consumption_total[!is.finite(d$ln_hh_consumption_total)] <- NA

# climate PCA (correlation matrix, like Stata pca) — span is what matters
pc_ok <- complete.cases(d[, c("TMP", "PER", "SNOW", "DEM", "lat")])
pc <- prcomp(d[pc_ok, c("TMP", "PER", "SNOW", "DEM", "lat")], center = TRUE, scale. = TRUE)
d[, paste0("climate_pc", 1:5)] <- NA
d[pc_ok, paste0("climate_pc", 1:5)] <- pc$x[, 1:5]

ctrl <- paste0("age + factor(sex) + factor(education) + ln_hh_consumption_total + ",
  "factor(retirement) + factor(pension_insurance) + factor(smoking_current) + ",
  "factor(alcohol_current) + climate_pc1 + climate_pc2 + climate_pc3 + climate_pc4 + ",
  "climate_pc5 + bmi + factor(hypertension) + factor(dyslipidemia) + factor(year)")

fit_one <- function(gvar, buf) {
  f <- as.formula(paste0("BAA_PhenoAge ~ factor(", gvar, ") + ", ctrl, " + (1 | id_num)"))
  m <- lmer(f, data = d, REML = FALSE)
  co <- summary(m)$coefficients
  rows <- list()
  for (lv in 2:4) {
    nm <- paste0("factor(", gvar, ")", lv)
    est <- co[nm, "Estimate"]; se <- co[nm, "Std. Error"]
    z <- est / se
    rows[[lv - 1]] <- data.frame(
      classification = "Original K-means (R replication)", buffer = buf,
      contrast = c("", "Low-EVI/High-O3 vs Reference", "High-EVI/Low-O3 vs Reference",
                   "High-EVI/High-O3 vs Reference")[lv],
      beta = round(est, 4), se = round(se, 4),
      pval = round(2 * pnorm(-abs(z)), 4),
      ci_low = round(est - 1.96 * se, 4), ci_high = round(est + 1.96 * se, 4),
      nobs = nobs(m))
  }
  do.call(rbind, rows)
}

out <- rbind(fit_one("joint_exp_1KM", "1KM"), fit_one("joint_exp_500m", "500m"))
write.csv(out, "output/tables/s20_unweighted_verification.csv", row.names = FALSE)
print(out)

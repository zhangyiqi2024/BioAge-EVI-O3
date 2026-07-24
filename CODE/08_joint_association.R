# =============================================================================
# 08_joint_association.R
# Primary joint-association analysis. Weighted linear mixed model
# (participant random intercept) via WeMix::mix
# with probability weights (cWeights = FALSE, matching Stata pweight) and
# robust (sandwich/Taylor) variance, matching Stata vce(robust).
# Adjusted marginal means and lincom-equivalent contrasts (including the
# difference-in-differences interaction) are computed by marginal
# standardization with the robust covariance matrix.
# Writes the official joint-association CSV outputs (joint_association_*_v4.csv).
# This script does not estimate or claim ozone offsetting/attenuation by EVI.
# =============================================================================

library(haven)
library(dplyr)
library(WeMix)

set.seed(20240101)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/08_joint_association.log", split = TRUE)
on.exit(sink(), add = TRUE)

df <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
df <- df[df$age > 45 & df$age < 90, ]
df$ln_hh_consumption_total <- ifelse(df$hh_consumption_total > 0,
                                     log(df$hh_consumption_total), NA_real_)

# Climate PCs (Stata: pca TMP PER SNOW DEM lat on the same filtered sample)
climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")
cc <- complete.cases(df[, climate_vars])
pc_fit <- prcomp(df[cc, climate_vars], center = TRUE, scale. = TRUE)
pcs <- predict(pc_fit, newdata = df[, climate_vars])
df[, paste0("climate_pc", 1:5)] <- pcs[, 1:5]

group_levels <- c("Low EVI / low O3", "Low EVI / high O3",
                  "High EVI / low O3", "High EVI / high O3")

results_list <- list()
means_list <- list()

for (s in c("1KM", "500m")) {
  joint <- paste0("joint_exp_", s)
  wt <- paste0("ipw_joint_", s)

  mvars <- c("BAA_PhenoAge", joint, wt, "id_num", "age", "sex", "education",
             "ln_hh_consumption_total", "retirement", "pension_insurance",
             "smoking_current", "alcohol_current", paste0("climate_pc", 1:5),
             "bmi", "hypertension", "dyslipidemia", "year")
  md <- df[, mvars]
  md <- md[complete.cases(md), ]
  md$joint <- factor(md[[joint]], levels = 1:4, labels = group_levels)
  md$id_num <- factor(md$id_num)
  md$year <- factor(md$year)
  for (v in c("sex", "education", "retirement", "pension_insurance",
              "smoking_current", "alcohol_current", "hypertension", "dyslipidemia"))
    md[[v]] <- factor(md[[v]])
  md$w <- md[[wt]]
  md$w2 <- 1  # unit weight at the participant (level-2) level

  f_fixed <- BAA_PhenoAge ~ joint + age + sex + education +
    ln_hh_consumption_total + retirement + pension_insurance +
    smoking_current + alcohol_current +
    climate_pc1 + climate_pc2 + climate_pc3 + climate_pc4 + climate_pc5 +
    bmi + hypertension + dyslipidemia + year

  fit <- mix(update(f_fixed, . ~ . + (1 | id_num)), data = md,
             weights = c("w", "w2"), cWeights = FALSE, verbose = FALSE)
  model_n <- nrow(md)
  cat("\n=== Scale:", s, "| N =", model_n, "===\n")

  beta <- coef(fit)                 # fixed-effect estimates
  V <- fit$cov_mat                  # robust (sandwich/Taylor) covariance
  beta <- beta[rownames(V)]         # align order

  # Marginal standardization: average model matrix per group level
  tt <- delete.response(terms(f_fixed))
  a_mat <- matrix(NA, nrow = 4, ncol = length(beta),
                  dimnames = list(group_levels, names(beta)))
  for (g in seq_along(group_levels)) {
    md_g <- md
    md_g$joint <- factor(group_levels[g], levels = group_levels)
    Xg <- model.matrix(tt, md_g)
    a_mat[g, ] <- colMeans(Xg)[names(beta)]
  }
  means <- as.numeric(a_mat %*% beta)
  se_means <- sqrt(diag(a_mat %*% V %*% t(a_mat)))

  for (j in 1:4) {
    means_list[[length(means_list) + 1]] <- data.frame(
      scale = s, group_order = j, exposure_group = group_levels[j],
      adjusted_mean = means[j], se = se_means[j],
      lci = means[j] - 1.96 * se_means[j], uci = means[j] + 1.96 * se_means[j],
      n = model_n)
  }

  ctr <- list(
    "Low EVI/high O3 vs reference"          = c(-1, 1, 0, 0),
    "High EVI/low O3 vs reference"          = c(-1, 0, 1, 0),
    "High EVI/high O3 vs reference"         = c(-1, 0, 0, 1),
    "O3 contrast within low EVI"            = c(-1, 1, 0, 0),
    "O3 contrast within high EVI"           = c(0, 0, -1, 1),
    # DiD = (m4 - m3) - (m2 - m1) = m1 - m2 - m3 + m4; the reference cell
    # (group 1, low EVI / low O3) enters with weight +1.
    "Difference-in-differences interaction" = c(1, -1, -1, 1)
  )
  for (nm in names(ctr)) {
    d <- as.numeric(t(ctr[[nm]]) %*% a_mat)   # contrast vector in beta space
    est <- sum(d * beta)
    se <- sqrt(as.numeric(t(d) %*% V %*% d))
    z <- est / se
    p <- 2 * pnorm(-abs(z))
    results_list[[length(results_list) + 1]] <- data.frame(
      scale = s, contrast = nm, estimate = est, se = se,
      lci = est - 1.96 * se, uci = est + 1.96 * se, p = p, n = model_n)
    cat(sprintf("%-42s est=%8.4f  se=%7.4f  p=%.4f\n", nm, est, se, p))
  }
}

contrasts_out <- do.call(rbind, results_list)
means_out <- do.call(rbind, means_list)

write.csv(contrasts_out, "output/tables/joint_association_contrasts_v4.csv",
          row.names = FALSE)
write.csv(means_out, "output/tables/joint_association_adjusted_means_v4.csv",
          row.names = FALSE)

cat("\n=== Adjusted means summary ===\n")
print(means_out)
cat("\nDone. CSVs written to output/tables/ (same filenames as Stata script).\n")

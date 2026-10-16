# =============================================================================
# 48_weight_diagnostics.R  (v5 revision, task T10b)
# IPW diagnostics for the joint-exposure weights (ipw_joint_1KM / ipw_joint_500m)
# on the authoritative joint-model frame (N=6,471):
#   - weight distribution (mean/sd/min/max/quantiles/CV)
#   - effective sample size (ESS = (sum w)^2 / sum w^2), overall and per group
#   - covariate balance: max pairwise SMD across the 4 joint groups,
#     unweighted vs weighted
# Inputs: DATA/analysis_ready/analysis_dataset.dta
# =============================================================================

suppressPackageStartupMessages({
  library(haven); library(dplyr)
})

set.seed(20240101)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/48_weight_diagnostics.log", split = TRUE)
on.exit(sink(), add = TRUE)

df <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
df <- df[df$age > 45 & df$age < 90, ]

# Joint-model frame mirrors 08_joint_association.R: complete on the joint
# grouping + weights + model covariates
df$ln_hh_consumption_total <- ifelse(df$hh_consumption_total > 0,
                                     log(df$hh_consumption_total), NA_real_)
cv <- c("TMP", "PER", "SNOW", "DEM", "lat")
cc <- complete.cases(df[, cv])
pc <- prcomp(df[cc, cv], center = TRUE, scale. = TRUE)
df[, paste0("climate_pc", 1:5)] <- predict(pc, newdata = df[, cv])[, 1:5]

covs <- c("age", "sex", "education", "ln_hh_consumption_total", "retirement",
          "pension_insurance", "smoking_current", "alcohol_current",
          paste0("climate_pc", 1:5), "bmi", "hypertension", "dyslipidemia", "year")

run_scale <- function(d, group_var, w_var, scale_lab) {
  d <- d[complete.cases(d[, c(group_var, w_var, covs, "BAA_PhenoAge")]), ]
  g <- factor(as.numeric(d[[group_var]]))   # drop haven_labelled class
  w <- d[[w_var]]
  cat("\n=== ", scale_lab, " | N =", nrow(d), " | groups:",
      paste(names(table(g)), table(g), sep = "=", collapse = ", "), "\n")

  qs <- quantile(w, c(.01, .05, .25, .5, .75, .95, .99))
  ess_all <- sum(w)^2 / sum(w^2)
  ess_g <- tapply(w, g, function(x) sum(x)^2 / sum(x^2))
  dist <- data.frame(
    scale = scale_lab, n = nrow(d), w_mean = mean(w), w_sd = sd(w),
    w_min = min(w), w_p01 = qs[1], w_p05 = qs[2], w_p25 = qs[3],
    w_median = qs[4], w_p75 = qs[5], w_p95 = qs[6], w_p99 = qs[7],
    w_max = max(w), w_cv = sd(w) / mean(w),
    ess = ess_all, ess_pct = 100 * ess_all / nrow(d),
    ess_g1 = ess_g[1], ess_g2 = ess_g[2], ess_g3 = ess_g[3], ess_g4 = ess_g[4])
  print(dist, digits = 4)

  # Max pairwise SMD per covariate, unweighted vs weighted
  gl <- names(table(g))
  bal <- lapply(covs, function(v) {
    x <- d[[v]]
    if (length(unique(x[!is.na(x)])) == 2) x <- as.numeric(as.character(x))
    x <- as.numeric(x)
    smd_pair <- function(wt) {
      prs <- combn(gl, 2, simplify = FALSE)
      mx <- 0
      for (pr in prs) {
        i1 <- g == pr[1]; i2 <- g == pr[2]
        if (is.null(wt)) {
          m1 <- mean(x[i1]); m2 <- mean(x[i2])
          v1 <- var(x[i1]);  v2 <- var(x[i2])
        } else {
          m1 <- weighted.mean(x[i1], wt[i1]); m2 <- weighted.mean(x[i2], wt[i2])
          v1 <- sum(wt[i1] * (x[i1] - m1)^2) / (sum(wt[i1]) - 1)
          v2 <- sum(wt[i2] * (x[i2] - m2)^2) / (sum(wt[i2]) - 1)
        }
        sp <- sqrt((v1 + v2) / 2)
        if (sp > 0) mx <- max(mx, abs(m1 - m2) / sp)
      }
      mx
    }
    data.frame(scale = scale_lab, covariate = v,
               max_smd_unweighted = smd_pair(NULL),
               max_smd_weighted = smd_pair(w))
  })
  bal <- do.call(rbind, bal)
  list(dist = dist, bal = bal)
}

r1 <- run_scale(df, "joint_exp_1KM", "ipw_joint_1KM", "1KM")
r2 <- run_scale(df, "joint_exp_500m", "ipw_joint_500m", "500m")

dist_all <- rbind(r1$dist, r2$dist)
bal_all <- rbind(r1$bal, r2$bal)
bal_all$imbalance_flag <- bal_all$max_smd_weighted > 0.1
cat("\nWeighted SMD > 0.1 flags:", sum(bal_all$imbalance_flag), "of", nrow(bal_all), "\n")
print(bal_all[order(-bal_all$max_smd_weighted), ][1:10, ], digits = 4)

write.csv(dist_all, "output/tables/weight_distribution_v5.csv", row.names = FALSE)
write.csv(bal_all, "output/tables/weight_balance_smd_v5.csv", row.names = FALSE)
cat("\nDone. output/tables/weight_distribution_v5.csv + weight_balance_smd_v5.csv\n")

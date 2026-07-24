# =============================================================================
# 30_iv_diagnostics.R
# R replication of 28_iv_analysis_stata.do / 29_iv_diagnostics_stata.do
# (Stata batch unavailable in this environment).
#
# Corrects three defects identified in the v4 audit:
#   (1) S11 endogeneity P-values were hardcoded / used an invalid variance
#       formula -> replaced with REAL control-function (Wooldridge robust
#       Durbin-Wu-Hausman) tests, cluster-robust by id_num.
#   (2) Dual-exposure first-stage F-statistics were copied from the
#       single-exposure models -> replaced with REAL Sanderson-Windmeijer
#       conditional F-statistics (cluster-robust Wald form) per endogenous
#       variable.
#   (3) Single-exposure F-statistics were hardcoded -> recomputed from the
#       actual first-stage regressions.
#
# Outputs:
#   output/tables/iv_validity_tests_v4.csv  (replaces S11 table body)
#   output/tables/iv_estimates_v4.csv       (2SLS + OLS coefficients)
#   logs/30_iv_diagnostics.log
# =============================================================================

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(fixest)
})

set.seed(20240101)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/30_iv_diagnostics.log", split = TRUE)
on.exit(sink(), add = TRUE)

# ---- data --------------------------------------------------------------------
df <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
df <- df[df$age > 45 & df$age < 90, ]
cat("Age-filtered sample N:", nrow(df), "\n")

df$ln_hh_consumption_total <- ifelse(df$hh_consumption_total > 0,
                                     log(df$hh_consumption_total), NA_real_)
df$log_EVI_1KM  <- log(df$EVI_1KM_mean4y + 0.001)
df$log_EVI_500m <- log(df$EVI_500m_mean4y + 0.001)
df$log_O3_1KM   <- log(df$O3_1KM_mean4y + 0.001)
df$log_O3_500m  <- log(df$O3_500m_mean4y + 0.001)
df$steep_slope_sq <- df$steep_slope_pct_25deg_500m_lag5^2

# Climate PCs (Stata: pca TMP PER SNOW DEM lat -> predict, score)
climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")
cc <- complete.cases(df[, climate_vars])
pc_fit <- prcomp(df[cc, climate_vars], center = TRUE, scale. = TRUE)
pcs <- predict(pc_fit, newdata = df[, climate_vars])
df[, paste0("climate_pc", 1:5)] <- pcs[, 1:5]
cat("PCA complete cases:", sum(cc), "\n")

# Controls (matches Stata spec exactly)
ctrl <- c("age", "factor(sex)", "factor(education)", "ln_hh_consumption_total",
          "factor(retirement)", "factor(pension_insurance)",
          "factor(smoking_current)", "factor(alcohol_current)",
          paste0("climate_pc", 1:5),
          "bmi", "factor(hypertension)", "factor(dyslipidemia)", "factor(year)")
ctrl_rhs <- paste(ctrl, collapse = " + ")

# Model frame shared by all IV models (complete cases on union of vars)
mvars <- c("BAA_PhenoAge", "id_num", "log_EVI_1KM", "log_EVI_500m",
           "log_O3_1KM", "log_O3_500m", "steep_slope_sq",
           "inversion_days1_lag5",
           "age", "sex", "education", "ln_hh_consumption_total", "retirement",
           "pension_insurance", "smoking_current", "alcohol_current",
           paste0("climate_pc", 1:5), "bmi", "hypertension", "dyslipidemia",
           "year")
md <- df[, mvars]
md <- md[complete.cases(md), ]
md$id_num <- factor(md$id_num)
cat("IV model frame N:", nrow(md), "\n\n")

# Helper: cluster-robust joint Wald F-test for named coefficients -------------
wald_F <- function(fit, coefs) {
  V <- vcov(fit)
  b <- coef(fit)[coefs]
  R <- diag(length(coefs))
  q <- length(coefs)
  stat <- as.numeric(t(R %*% b) %*% solve(R %*% V[coefs, coefs] %*% t(R)) %*%
                       (R %*% b)) / q
  list(F = stat, p = pf(stat, q, df.residual(fit), lower.tail = FALSE))
}

# =============================================================================
# Part 1: single-exposure models (4) — first-stage F, OLS, 2SLS, DWH
# =============================================================================
single_specs <- list(
  list(exp = "log_EVI_1KM",  z = "steep_slope_sq",       label = "EVI 1KM"),
  list(exp = "log_EVI_500m", z = "steep_slope_sq",       label = "EVI 500m"),
  list(exp = "log_O3_1KM",   z = "inversion_days1_lag5", label = "O3 1KM"),
  list(exp = "log_O3_500m",  z = "inversion_days1_lag5", label = "O3 500m")
)

rows_validity <- list()
rows_est <- list()

for (sp in single_specs) {
  cat("=== Single-exposure:", sp$label, "===\n")

  # First stage
  f_fs <- as.formula(paste(sp$exp, "~", sp$z, "+", ctrl_rhs))
  fs <- feols(f_fs, data = md, cluster = ~id_num)
  wt <- wald_F(fs, sp$z)
  cat(sprintf("  First-stage F = %.2f (p = %.4g)\n", wt$F, wt$p))

  # OLS
  f_ols <- as.formula(paste("BAA_PhenoAge ~", sp$exp, "+", ctrl_rhs))
  ols <- feols(f_ols, data = md, cluster = ~id_num)

  # 2SLS
  f_iv <- as.formula(paste0("BAA_PhenoAge ~ ", ctrl_rhs, " | ",
                            sp$exp, " ~ ", sp$z))
  iv <- feols(f_iv, data = md, cluster = ~id_num)

  # Control-function (Wooldridge robust DWH) endogeneity test
  md$resid_fs <- residuals(fs)
  f_cf <- as.formula(paste("BAA_PhenoAge ~", sp$exp, "+ resid_fs +", ctrl_rhs))
  cf <- feols(f_cf, data = md, cluster = ~id_num)
  dwh_p <- pvalue(cf)["resid_fs"]
  cat(sprintf("  DWH (control-function) p = %.4g\n", dwh_p))
  md$resid_fs <- NULL

  rows_validity[[length(rows_validity) + 1]] <- data.frame(
    exposure = ifelse(grepl("EVI", sp$exp), "log(EVI)", "log(O3)"),
    model = "Single", buffer = ifelse(grepl("1KM", sp$exp), "1KM", "500m"),
    instrument = ifelse(grepl("EVI", sp$exp), "Steep slope squared",
                        "Temp. inversion freq."),
    F_stat = round(wt$F, 2), F_p = signif(wt$p, 3),
    weak_iv = ifelse(wt$F > 10, "Strong", ifelse(wt$F > 5, "Moderate", "Weak")),
    endog_p = signif(dwh_p, 3),
    endog_conclusion = ifelse(dwh_p < 0.10, "Present", "Not detected"),
    identification = "Exact")

  ci <- confint(iv)[paste0("fit_", sp$exp), ]
  rows_est[[length(rows_est) + 1]] <- data.frame(
    exposure = ifelse(grepl("EVI", sp$exp), "EVI", "O3"),
    model = "Single-exposure",
    buffer = ifelse(grepl("1KM", sp$exp), "1KM", "500m"),
    coef_iv = coef(iv)[paste0("fit_", sp$exp)],
    se_iv = se(iv)[paste0("fit_", sp$exp)],
    p_iv = pvalue(iv)[paste0("fit_", sp$exp)],
    lci_iv = ci[1], uci_iv = ci[2],
    coef_ols = coef(ols)[sp$exp], se_ols = se(ols)[sp$exp],
    n = nobs(iv))
}

# =============================================================================
# Part 2: dual-exposure models (2) — SW conditional F, 2SLS, joint DWH
# =============================================================================
# Sanderson-Windmeijer conditional F (cluster-robust Wald form):
#   FWL: partial out exogenous controls W from x1, x2, z1, z2;
#   fitted values of the OTHER endogenous regressor on both instruments;
#   partial those fitted values out of x1 and both instruments;
#   cluster-robust joint Wald that the remaining instrument variation = 0.
sw_cond_F <- function(data, x1, x2, z1, z2) {
  f_partial <- function(v) {
    residuals(feols(as.formula(paste(v, "~", ctrl_rhs)), data = data))
  }
  x1s <- f_partial(x1); x2s <- f_partial(x2)
  z1s <- f_partial(z1); z2s <- f_partial(z2)
  x2hat <- fitted(lm(x2s ~ z1s + z2s))
  resid_on <- function(y) residuals(lm(y ~ x2hat))
  d <- data.frame(x1t = resid_on(x1s), z1t = resid_on(z1s),
                  z2t = resid_on(z2s), cl = data$id_num)
  fit <- feols(x1t ~ z1t + z2t - 1, data = d, cluster = ~cl)
  kept <- intersect(c("z1t", "z2t"), names(coef(fit)))
  wald_F(fit, kept)
}

dual_specs <- list(
  list(evi = "log_EVI_1KM",  o3 = "log_O3_1KM",  buffer = "1KM"),
  list(evi = "log_EVI_500m", o3 = "log_O3_500m", buffer = "500m")
)

for (sp in dual_specs) {
  cat("\n=== Dual-exposure:", sp$buffer, "===\n")

  # 2SLS with two endogenous regressors
  f_iv <- as.formula(paste0("BAA_PhenoAge ~ ", ctrl_rhs, " | ",
                            sp$evi, " + ", sp$o3,
                            " ~ steep_slope_sq + inversion_days1_lag5"))
  iv <- feols(f_iv, data = md, cluster = ~id_num)

  # Sanderson-Windmeijer conditional F per endogenous variable
  sw_evi <- sw_cond_F(md, sp$evi, sp$o3, "steep_slope_sq", "inversion_days1_lag5")
  sw_o3  <- sw_cond_F(md, sp$o3, sp$evi, "inversion_days1_lag5", "steep_slope_sq")
  cat(sprintf("  SW conditional F: EVI = %.2f (p = %.4g), O3 = %.2f (p = %.4g)\n",
              sw_evi$F, sw_evi$p, sw_o3$F, sw_o3$p))

  # Control-function DWH: two first stages, both residuals in outcome equation
  f_fs1 <- as.formula(paste(sp$evi, "~ steep_slope_sq + inversion_days1_lag5 +",
                            ctrl_rhs))
  f_fs2 <- as.formula(paste(sp$o3, "~ steep_slope_sq + inversion_days1_lag5 +",
                            ctrl_rhs))
  md$r1 <- residuals(feols(f_fs1, data = md, cluster = ~id_num))
  md$r2 <- residuals(feols(f_fs2, data = md, cluster = ~id_num))
  f_cf <- as.formula(paste("BAA_PhenoAge ~", sp$evi, "+", sp$o3,
                           "+ r1 + r2 +", ctrl_rhs))
  cf <- feols(f_cf, data = md, cluster = ~id_num)
  dwh_p_evi <- pvalue(cf)["r1"]
  dwh_p_o3  <- pvalue(cf)["r2"]
  dwh_joint <- wald_F(cf, c("r1", "r2"))
  cat(sprintf("  DWH p: EVI resid = %.4g, O3 resid = %.4g, joint F = %.3f (p = %.4g)\n",
              dwh_p_evi, dwh_p_o3, dwh_joint$F, dwh_joint$p))

  for (which in c("EVI", "O3")) {
    x <- if (which == "EVI") sp$evi else sp$o3
    sw <- if (which == "EVI") sw_evi else sw_o3
    dp <- if (which == "EVI") dwh_p_evi else dwh_p_o3
    rows_validity[[length(rows_validity) + 1]] <- data.frame(
      exposure = paste0("log(", which, ")"),
      model = "Dual", buffer = sp$buffer,
      instrument = "Steep slope + Temp. inversion",
      F_stat = round(sw$F, 2), F_p = signif(sw$p, 3),
      weak_iv = ifelse(sw$F > 10, "Strong",
                       ifelse(sw$F > 5, "Moderate", "Weak")),
      endog_p = signif(dp, 3),
      endog_conclusion = ifelse(dp < 0.10, "Present", "Not detected"),
      identification = "Exact")

    ci <- confint(iv)[paste0("fit_", x), ]
    rows_est[[length(rows_est) + 1]] <- data.frame(
      exposure = which, model = "Dual-exposure", buffer = sp$buffer,
      coef_iv = coef(iv)[paste0("fit_", x)],
      se_iv = se(iv)[paste0("fit_", x)],
      p_iv = pvalue(iv)[paste0("fit_", x)],
      lci_iv = ci[1], uci_iv = ci[2],
      coef_ols = NA_real_, se_ols = NA_real_,
      n = nobs(iv))
  }
  md$r1 <- NULL; md$r2 <- NULL
}

# ---- write outputs -----------------------------------------------------------
validity <- do.call(rbind, rows_validity)
ests <- do.call(rbind, rows_est)

write.csv(validity, "output/tables/iv_validity_tests_v4.csv", row.names = FALSE)
write.csv(ests, "output/tables/iv_estimates_v4.csv", row.names = FALSE)

cat("\n=== IV validity tests (v4, recomputed) ===\n")
print(validity, row.names = FALSE)
cat("\n=== 2SLS estimates (v4, recomputed) ===\n")
print(ests, row.names = FALSE)

cat("\nNotes for SI:\n")
cat("- Single-exposure F: cluster-robust F-test of the instrument in the first\n")
cat("  stage (equals squared t-stat for a single instrument).\n")
cat("- Dual-exposure F: Sanderson-Windmeijer conditional F (cluster-robust Wald\n")
cat("  form), testing instrument relevance conditional on the other endogenous\n")
cat("  regressor.\n")
cat("- Endogeneity: control-function (Wooldridge robust Durbin-Wu-Hausman)\n")
cat("  test; p-value of the first-stage residual(s) in the outcome equation,\n")
cat("  cluster-robust by participant.\n")
cat("- All models exactly identified; overidentification tests not applicable.\n")
cat("\nDone. CSVs written to output/tables/.\n")

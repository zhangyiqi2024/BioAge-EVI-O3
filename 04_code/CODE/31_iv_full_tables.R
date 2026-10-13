# =============================================================================
# 31_iv_full_tables.R
# Full-coefficient export of the six v4 IV models (same specs/frame as
# 30_iv_diagnostics.R) for regenerating the 06_tables presentation files,
# plus a forest plot of the eight exposure coefficients.
#
# Outputs:
#   output/tables/iv_full_coefficients_v4.csv  (long: column, variable, beta, se, p)
#   output/figures/iv_forest_v4.png            (2000x1455, replaces stale embedded figure)
#   logs/31_iv_full_tables.log
# =============================================================================

suppressPackageStartupMessages({ library(haven); library(dplyr); library(fixest) })

set.seed(20240101)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/31_iv_full_tables.log", split = TRUE)
on.exit(sink(), add = TRUE)

# ---- data (identical to 10b) --------------------------------------------------
df <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
df <- df[df$age > 45 & df$age < 90, ]
df$ln_hh_consumption_total <- ifelse(df$hh_consumption_total > 0,
                                     log(df$hh_consumption_total), NA_real_)
df$log_EVI_1KM  <- log(df$EVI_1KM_mean4y + 0.001)
df$log_EVI_500m <- log(df$EVI_500m_mean4y + 0.001)
df$log_O3_1KM   <- log(df$O3_1KM_mean4y + 0.001)
df$log_O3_500m  <- log(df$O3_500m_mean4y + 0.001)
df$steep_slope_sq <- df$steep_slope_pct_25deg_500m_lag5^2

climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")
cc <- complete.cases(df[, climate_vars])
pc_fit <- prcomp(df[cc, climate_vars], center = TRUE, scale. = TRUE)
df[, paste0("climate_pc", 1:5)] <- predict(pc_fit, newdata = df[, climate_vars])[, 1:5]

ctrl <- c("age", "factor(sex)", "factor(education)", "ln_hh_consumption_total",
          "factor(retirement)", "factor(pension_insurance)",
          "factor(smoking_current)", "factor(alcohol_current)",
          paste0("climate_pc", 1:5),
          "bmi", "factor(hypertension)", "factor(dyslipidemia)", "factor(year)")
ctrl_rhs <- paste(ctrl, collapse = " + ")

mvars <- c("BAA_PhenoAge", "id_num", "log_EVI_1KM", "log_EVI_500m",
           "log_O3_1KM", "log_O3_500m", "steep_slope_sq", "inversion_days1_lag5",
           "age", "sex", "education", "ln_hh_consumption_total", "retirement",
           "pension_insurance", "smoking_current", "alcohol_current",
           paste0("climate_pc", 1:5), "bmi", "hypertension", "dyslipidemia", "year")
md <- df[complete.cases(df[, mvars]), mvars]
md$id_num <- factor(md$id_num)
cat("IV model frame N:", nrow(md), "\n")

# ---- fit the six models -------------------------------------------------------
fit_iv <- function(rhs_exo) {
  feols(as.formula(paste0("BAA_PhenoAge ~ ", ctrl_rhs, " | ", rhs_exo)),
        data = md, cluster = ~id_num)
}
fits <- list(
  `EVI_1KM`  = fit_iv("log_EVI_1KM ~ steep_slope_sq"),
  `EVI_500m` = fit_iv("log_EVI_500m ~ steep_slope_sq"),
  `O3_1KM`   = fit_iv("log_O3_1KM ~ inversion_days1_lag5"),
  `O3_500m`  = fit_iv("log_O3_500m ~ inversion_days1_lag5"),
  `Dual_1KM` = fit_iv("log_EVI_1KM + log_O3_1KM ~ steep_slope_sq + inversion_days1_lag5"),
  `Dual_500m`= fit_iv("log_EVI_500m + log_O3_500m ~ steep_slope_sq + inversion_days1_lag5")
)

# ---- long full-coefficient table ----------------------------------------------
rows <- list()
for (nm in names(fits)) {
  f <- fits[[nm]]
  V <- summary(f)$coeftable
  for (v in rownames(V)) {
    rows[[length(rows) + 1]] <- data.frame(
      column = nm, variable = v,
      beta = V[v, "Estimate"], se = V[v, "Std. Error"], p = V[v, "Pr(>|t|)"],
      n = nobs(f))
  }
}
full <- do.call(rbind, rows)
write.csv(full, "output/tables/iv_full_coefficients_v4.csv", row.names = FALSE)
cat("full coefficient rows:", nrow(full), "\n")

# ---- verify exposure rows against iv_estimates_v4.csv -------------------------
ref <- read.csv("output/tables/iv_estimates_v4.csv")
lookup <- full[grepl("^fit_log", full$variable), ]
lookup$var <- gsub("^fit_log_", "", lookup$variable)
lookup$key <- paste(lookup$column, lookup$var, sep = "||")
ref$var <- paste0(ref$exposure, "_", ref$buffer)
ref$col <- ifelse(ref$model == "Dual-exposure",
                  paste0("Dual_", ref$buffer), paste0(ref$exposure, "_", ref$buffer))
ref$key <- paste(ref$col, ref$var, sep = "||")
mrg <- merge(lookup, ref, by = "key")
stopifnot(nrow(mrg) == 8)
tol_beta <- max(abs(mrg$beta - mrg$coef_iv))
tol_se   <- max(abs(mrg$se - mrg$se_iv))
cat(sprintf("verification: max |dBeta| = %.3g, max |dSE| = %.3g (must be < 1e-6 / 1e-6)\n",
            tol_beta, tol_se))
stopifnot(tol_beta < 1e-6, tol_se < 1e-6)
cat("VERIFIED: full tables consistent with iv_estimates_v4.csv\n")

# ---- forest plot from authoritative v4 estimates ------------------------------
est <- read.csv("output/tables/iv_estimates_v4.csv")
est <- est[order(match(est$model, c("Single-exposure", "Dual-exposure")),
                 match(est$exposure, c("EVI", "O3")), match(est$buffer, c("1KM", "500m"))), ]
est$label <- paste(est$exposure, ifelse(est$model == "Single-exposure", "Single", "Dual"), est$buffer)
ylab <- rev(est$label)

png("output/figures/iv_forest_v4.png", width = 2000, height = 1455, res = 300)
par(mar = c(4.5, 7.5, 2.5, 1.5), family = "sans")
y <- rev(seq_len(nrow(est)))
cols <- ifelse(est$exposure == "EVI", "#1a9850", "#d73027")
pchs <- ifelse(est$model == "Single-exposure", 16, 17)
plot(NA, xlim = range(c(est$X2.5.., est$X97.5..)) * 1.04, ylim = c(0.4, nrow(est) + 0.6),
     yaxt = "n", xlab = "2SLS coefficient (years per log-unit)", ylab = "",
     main = "Forest Plot of Instrumental Variable Estimates", cex.main = 1.1)
abline(v = 0, lty = 2, col = "grey50")
segments(est$X2.5.., y, est$X97.5.., y, col = cols, lwd = 2)
points(est$coef_iv, y, pch = pchs, col = cols, cex = 1.3)
axis(2, at = y, labels = est$label, las = 1, cex.axis = 0.85)
legend("bottomleft", legend = c("EVI (single)", "EVI (dual)", "O3 (single)", "O3 (dual)"),
       col = c("#1a9850", "#1a9850", "#d73027", "#d73027"), pch = c(16, 17, 16, 17),
       bty = "n", cex = 0.85, inset = 0.02)
dev.off()
cat("forest plot saved: output/figures/iv_forest_v4.png\n")

# =============================================================================
# 40_clinical_benchmarks.R
# Recomputes Supplementary Table S21 (clinical benchmarking of environmental
# vs behavioural BAA associations) for the v4 corrective revision.
#
# Corrections vs the old S21:
#   (1) Environmental rows previously used the pooled meta-analytic estimates
#       (-1.283 / +1.000), which were REMOVED in v4 -> replaced with the
#       mutually adjusted (dual-exposure) weighted mixed-model estimates
#       from Exposure_Effects_Summary.csv (EVI -1.318; O3 +1.114, 1 km).
#   (2) Behavioural rows previously disagreed between the SI text block
#       (smoking +0.83, P=0.120) and the S21 table (+0.442, P=0.269) ->
#       recomputed once, with the linear mixed-effects specification
#       described in the S21 methods text, and used in BOTH places.
#   (3) Sample counts updated to the current analytical dataset.
#
# Outputs:
#   output/tables/s21_clinical_benchmarks_v4.csv
#   logs/40_clinical_benchmarks.log
# =============================================================================

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(lme4)
})

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/40_clinical_benchmarks.log", split = TRUE)
on.exit(sink(), add = TRUE)

df <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
df <- df[df$age > 45 & df$age < 90, ]
df$ln_hh_consumption_total <- ifelse(df$hh_consumption_total > 0,
                                     log(df$hh_consumption_total), NA_real_)
df$obese   <- as.integer(df$bmi >= 28)
df$low_edu <- as.integer(df$education == 1)   # less than junior high school

climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")
cc <- complete.cases(df[, climate_vars])
pc_fit <- prcomp(df[cc, climate_vars], center = TRUE, scale. = TRUE)
df[, paste0("climate_pc", 1:5)] <- predict(pc_fit, newdata = df[, climate_vars])[, 1:5]

base_covs <- c("age", "sex", "education", "ln_hh_consumption_total",
               "retirement", "pension_insurance", "smoking_current",
               "alcohol_current", "bmi", "hypertension", "dyslipidemia",
               paste0("climate_pc", 1:5), "year")

mvars <- c("BAA_PhenoAge", "id_num", "obese", "low_edu", base_covs)
md <- df[, mvars]
md <- md[complete.cases(md), ]
md$id_num <- factor(md$id_num)
for (v in c("sex", "education", "retirement", "pension_insurance",
            "smoking_current", "alcohol_current", "hypertension",
            "dyslipidemia", "year"))
  md[[v]] <- factor(md[[v]])
cat("Benchmark model frame: N obs =", nrow(md),
    " N participants =", length(unique(md$id_num)), "\n\n")

# Each risk factor enters alone; its own domain is dropped from adjustment.
specs <- list(
  list(var = "smoking_current", label = "Current smoking (vs. never)",
       drop = "smoking_current"),
  list(var = "obese",           label = "Obesity (BMI >= 28)",
       drop = "bmi"),
  list(var = "hypertension",    label = "Hypertension",
       drop = "hypertension"),
  list(var = "low_edu",         label = "Low education (< junior high)",
       drop = "education"),
  list(var = "alcohol_current", label = "Current alcohol use",
       drop = "alcohol_current")
)

rows <- list()
for (sp in specs) {
  covs <- setdiff(base_covs, sp$drop)
  f <- as.formula(paste("BAA_PhenoAge ~", sp$var, "+",
                        paste(covs, collapse = " + "), "+ (1 | id_num)"))
  fit <- lmer(f, data = md, REML = FALSE)
  cf <- summary(fit)$coefficients
  # coefficient row: for factor predictors it is the non-reference level
  rn <- rownames(cf)
  hit <- grep(paste0("^", sp$var), rn, value = TRUE)
  est <- cf[hit[1], "Estimate"]; se <- cf[hit[1], "Std. Error"]
  z <- est / se; p <- 2 * pnorm(-abs(z))
  rows[[length(rows) + 1]] <- data.frame(
    domain = "Behavioural", risk_factor = sp$label,
    baa_diff = round(est, 3), lci = round(est - 1.96 * se, 3),
    uci = round(est + 1.96 * se, 3), p = signif(p, 3),
    significant = ifelse(p < 0.05, "Yes", "No"))
  cat(sprintf("%-32s %+7.3f [%+7.3f, %+7.3f]  p = %.4g\n",
              sp$label, est, est - 1.96*se, est + 1.96*se, p))
}

# Environmental rows: mutually adjusted (dual-exposure) WME estimates, 1 km
env <- read.csv("output/tables/Exposure_Effects_Summary.csv")
env_rows <- list(
  list(id = "Joint_EVI_1KM_WME", label = "High vs Low EVI (dual-exposure, 1 km)"),
  list(id = "Joint_O3_1KM_WME",  label = "High vs Low O3 (dual-exposure, 1 km)"))
for (er in env_rows) {
  r <- env[env$Model == er$id, ]
  coef <- as.numeric(sub("\\*+", "", r$Coefficient))
  ci <- gsub("[()]", "", r$X95..CI); ci <- as.numeric(strsplit(ci, ",")[[1]])
  z <- coef / r$SE; p <- 2 * pnorm(-abs(z))
  rows[[length(rows) + 1]] <- data.frame(
    domain = "Environmental", risk_factor = er$label,
    baa_diff = round(coef, 3), lci = round(ci[1], 3), uci = round(ci[2], 3),
    p = signif(p, 3), significant = ifelse(p < 0.05, "Yes", "No"))
  cat(sprintf("%-32s %+7.3f [%+7.3f, %+7.3f]  p = %.4g\n",
              er$label, coef, ci[1], ci[2], p))
}

# IQR rows: current primary mutually adjusted continuous specification.
primary_file <- "output/tables/primary_continuous_models_v5.csv"
if(!file.exists(primary_file)) {
  rscript <- file.path(R.home("bin"), if(.Platform$OS.type=="windows") "Rscript.exe" else "Rscript")
  status <- system2(rscript, "CODE/43_primary_continuous_models.R")
  if(status!=0) stop("Unable to generate primary continuous model table.")
}
primary <- read.csv(primary_file)
primary <- primary[primary$scale=="1KM" & primary$weight_spec=="unweighted" &
                   primary$model=="dual-exposure",]
stopifnot(nrow(primary)==2)
for(i in seq_len(nrow(primary))) {
  r <- primary[i,]
  rows[[length(rows)+1]] <- data.frame(
    domain="Environmental (IQR)",
    risk_factor=paste("Per IQR",r$exposure,"(dual-exposure, 1 km)"),
    baa_diff=round(r$beta_iqr,3),lci=round(r$lci_iqr,3),uci=round(r$uci_iqr,3),
    p=signif(r$p,3),significant=ifelse(r$p<0.05,"Yes","No"))
}

out <- do.call(rbind, rows)
out$n_obs <- nrow(md)
out$n_participants <- length(unique(md$id_num))
iqr_rows <- which(out$domain=="Environmental (IQR)")
out$n_obs[iqr_rows] <- primary$n
out$n_participants[iqr_rows] <- primary$n_participants
write.csv(out, "output/tables/s21_clinical_benchmarks_v4.csv", row.names = FALSE)

cat("\n=== S21 v4 (recomputed) ===\n")
print(out, row.names = FALSE)
cat("\nDone. CSV written to output/tables/s21_clinical_benchmarks_v4.csv\n")

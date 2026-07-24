# =============================================================================
# 41_attrition_analysis.R
# v4 corrective revision — completes three analyses that were cited in the
# manuscript but not delivered:
#   (1) Supplementary Table S18: baseline (2011) characteristics of
#       participants excluded for prevalent major chronic disease vs the
#       retained analytical sample, with standardised differences.
#   (2) Attrition comparison: baseline characteristics of participants
#       retained in the analytical cohort vs eligible 2011 participants
#       who were not retained ("2011-only").
#   (3) IPAW: primary joint EVI x O3 models re-fitted with combined
#       inverse-probability-of-attrition x exposure weights (primary-scale
#       BAA outcome, corrected exposure weights).
#   (4) Pattern-mixture sensitivity: delta = 0,1,2,3 years. Because the
#       primary-analysis BAA derivation is not fully documented in the
#       delivered code (its linear predictor cannot be reproduced for
#       participants outside the analytical sample), this analysis is run
#       on the documented recalibrated PhenoAge scale (post-QC reference
#       file), which correlates ~0.68 with the primary BAA. All rows
#       (stayers + attriters) use that single scale; the delta = 0 run is
#       the internal reference. Attriters require the full model covariate
#       set, so this analysis uses covariate-complete attriters.
# Also recomputes the honest sequential cohort-construction flow (S4 Panel A).
#
# Data provenance (raw CHARLS derivatives, NOT redistributable; paths are
# relative to the project root, three levels above this script's run dir):
#   ../../../DATA/原始数据及样本筛选/data_生物年龄.dta          (screening input)
#   ../../../DATA/原始数据及样本筛选/data_PhenoAge_计算完成.dta  (post-QC, 13,111)
#   ../../../DATA/Data_LAST总.dta                              (master merge)
# Analytical sample: DATA/analysis_ready/analysis_dataset.dta (corrected IPW).
#
# Outputs (all under output/tables/):
#   s18_excluded_vs_included_v4.csv
#   attrition_baseline_comparison_v4.csv
#   attrition_ipaw_joint_v4.csv
#   attrition_pattern_mixture_v4.csv
#   s4_cohort_flow_v4.csv
# Log: logs/41_attrition_analysis.log
# =============================================================================

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(WeMix)
})

set.seed(20240101)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/41_attrition_analysis.log", split = TRUE)
on.exit(sink(), add = TRUE)

ROOT <- "../../.."   # project root relative to Submission/final_Submission/04_code

# ---- 0. analytical sample ----------------------------------------------------
ana <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
ana <- ana[ana$age > 45 & ana$age < 90, ]
ana11 <- ana[ana$year == 2011, ]
cat("Analytical:", nrow(ana), "obs,", length(unique(ana$id_num)),
    "participants; baseline rows:", nrow(ana11), "\n")
# ana$id is the CHARLS string id (verified identical to the id<->id_num bridge)

# ---- 1. honest S4 cohort flow -------------------------------------------------
cat("\n=== S4 cohort flow recompute ===\n")
raw <- as.data.frame(read_dta(file.path(
  ROOT, "DATA/原始数据及样本筛选/data_生物年龄.dta")))
n0 <- nrow(raw)
s1 <- raw[!is.na(raw$city), ]
s2 <- s1[s1$age >= 45, ]
s3 <- s2[s2$wave %in% c(1, 3), ]
biomarkers <- c("bl_wbc", "bl_mcv", "bl_bun", "bl_crea", "bl_crp",
                "bl_cho", "bl_hbalc", "bl_ua", "sbp", "bl_glu")
s4 <- s3[rowSums(is.na(s3[, biomarkers])) == 0, ]
diseases <- c("癌症", "糖尿病", "心脏病", "中风", "肺病", "记忆疾病")
w1_dis <- s4$wave == 1 & rowSums(s4[, diseases] == 1, na.rm = TRUE) > 0
dis_ids <- s4$id[w1_dis]
s5 <- s4[!s4$id %in% dis_ids, ]
ref <- as.data.frame(read_dta(file.path(
  ROOT, "DATA/原始数据及样本筛选/data_PhenoAge_计算完成.dta")))
both_ids <- intersect(ref$id[ref$wave == 1], ref$id[ref$wave == 3])
s7 <- ref[ref$id %in% both_ids, ]

flow <- data.frame(
  step = c("0 Original merged observations (CHARLS waves 1-5)",
           "1 Non-missing residence (city)",
           "2 Age >= 45 years",
           "3 Biomarker waves (2011 and 2015)",
           "4 Complete core biomarkers",
           "5 No prevalent major chronic disease at baseline",
           "6 PhenoAge quality control (fasting, ranges, CRP < 10)",
           "7 Observed at both biomarker waves",
           "8 Analytical sample (45 < age < 90, exposures + covariates)"),
  n_obs = c(n0, nrow(s1), nrow(s2), nrow(s3), nrow(s4), nrow(s5),
            nrow(ref), nrow(s7), nrow(ana)),
  n_participants = c(rep(NA_integer_, 7), length(unique(s7$id)),
                     length(unique(ana$id_num))))
print(flow, row.names = FALSE)
write.csv(flow, "output/tables/s4_cohort_flow_v4.csv", row.names = FALSE)

# ---- 2. master-merge covariates + exposures ----------------------------------
need_cols <- c("id", "year", "性别", "肥胖指数bmi", "是否吸烟", "是否饮酒",
               "高血压", "血脂异常", "教育统一分类", "农村或城市", "婚姻",
               "EVI_1KM", "O3_1KM", "EVI_500m", "O3_500m",
               "EVI_1KM_mean4y", "O3_1KM_mean4y",
               "EVI_500m_lag1", "EVI_500m_lag2", "EVI_500m_lag3",
               "EVI_500m_lag4", "O3_500m_lag1", "O3_500m_lag2",
               "O3_500m_lag3", "O3_500m_lag4",
               "TMP", "PER", "SNOW", "DEM", "lat",
               "家庭总消费", "是否退休", "养老保险")
master <- as.data.frame(read_dta(file.path(ROOT, "DATA/Data_LAST总.dta"),
                                 col_select = any_of(need_cols)))
rename_map <- c("性别" = "sex_c", "肥胖指数bmi" = "bmi",
                "是否吸烟" = "smoking_current",
                "是否饮酒" = "alcohol_current", "高血压" = "hypertension",
                "血脂异常" = "dyslipidemia", "教育统一分类" = "education",
                "农村或城市" = "urban_rural", "婚姻" = "marital_status",
                "家庭总消费" = "hh_consumption_total",
                "是否退休" = "retirement", "养老保险" = "pension_insurance")
for (cn in intersect(names(rename_map), names(master)))
  names(master)[names(master) == cn] <- rename_map[cn]

# master EVI is on a 0-10 scale; the analytical EVI is 0-1 -> rescale
for (v in c("EVI_1KM", "EVI_500m", "EVI_1KM_mean4y"))
  master[[v]] <- master[[v]] / 10
# 500-m mean4y (mean of the four preceding years, same window as the 1-km
# mean4y; verified cor = 0.99 with the delivered 1-km mean4y)
master$EVI_500m_mean4y <-
  rowMeans(master[, paste0("EVI_500m_lag", 1:4)], na.rm = FALSE) / 10
master$O3_500m_mean4y <-
  rowMeans(master[, paste0("O3_500m_lag", 1:4)], na.rm = FALSE)

m11 <- master[master$year == 2011, ]
cat("Master rows: 2011 =", nrow(m11), "\n")

# ---- helpers ------------------------------------------------------------------
stat_row <- function(x_e, x_i, binary = FALSE) {
  if (binary) {
    pe <- mean(x_e, na.rm = TRUE); pi_ <- mean(x_i, na.rm = TRUE)
    p <- (pe * sum(!is.na(x_e)) + pi_ * sum(!is.na(x_i))) /
         (sum(!is.na(x_e)) + sum(!is.na(x_i)))
    d <- (pe - pi_) / sqrt(p * (1 - p))
    list(e = sprintf("%d (%.1f)", round(pe * sum(!is.na(x_e))), 100 * pe),
         i = sprintf("%d (%.1f)", round(pi_ * sum(!is.na(x_i))), 100 * pi_),
         d = round(d, 3),
         n_e = sum(!is.na(x_e)), n_i = sum(!is.na(x_i)))
  } else {
    me <- mean(x_e, na.rm = TRUE); mi <- mean(x_i, na.rm = TRUE)
    sp <- sqrt((var(x_e, na.rm = TRUE) + var(x_i, na.rm = TRUE)) / 2)
    list(e = sprintf("%.2f (%.2f)", me, sd(x_e, na.rm = TRUE)),
         i = sprintf("%.2f (%.2f)", mi, sd(x_i, na.rm = TRUE)),
         d = round((me - mi) / sp, 3),
         n_e = sum(!is.na(x_e)), n_i = sum(!is.na(x_i)))
  }
}

group_levels <- c("Low EVI / low O3", "Low EVI / high O3",
                  "High EVI / low O3", "High EVI / high O3")
ctr <- list(
  "Low EVI/high O3 vs reference"          = c(-1, 1, 0, 0),
  "High EVI/low O3 vs reference"          = c(-1, 0, 1, 0),
  "High EVI/high O3 vs reference"         = c(-1, 0, 0, 1),
  "O3 contrast within high EVI"           = c(0, 0, -1, 1),
  # DiD = (m4 - m3) - (m2 - m1) = m1 - m2 - m3 + m4; the reference cell
  # (group 1, low EVI / low O3) enters with weight +1.
  "Difference-in-differences interaction" = c(1, -1, -1, 1))

fit_joint <- function(data, joint_var, wt_var, outcome = "BAA_PhenoAge") {
  mvars <- c(outcome, joint_var, wt_var, "id_num", "age", "sex",
             "education", "ln_hh_consumption_total", "retirement",
             "pension_insurance", "smoking_current", "alcohol_current",
             paste0("climate_pc", 1:5), "bmi", "hypertension",
             "dyslipidemia", "year")
  md <- data[complete.cases(data[, mvars]), mvars]
  md$joint <- factor(md[[joint_var]], levels = 1:4, labels = group_levels)
  md$y <- md[[outcome]]
  md$id_num <- factor(md$id_num); md$year <- factor(md$year)
  for (v in c("sex", "education", "retirement", "pension_insurance",
              "smoking_current", "alcohol_current", "hypertension",
              "dyslipidemia"))
    md[[v]] <- factor(md[[v]])
  md$w <- md[[wt_var]]; md$w2 <- 1
  f_fixed <- y ~ joint + age + sex + education +
    ln_hh_consumption_total + retirement + pension_insurance +
    smoking_current + alcohol_current + climate_pc1 + climate_pc2 +
    climate_pc3 + climate_pc4 + climate_pc5 + bmi + hypertension +
    dyslipidemia + year
  fit <- mix(update(f_fixed, . ~ . + (1 | id_num)), data = md,
             weights = c("w", "w2"), cWeights = FALSE, verbose = FALSE)
  beta <- coef(fit); V <- fit$cov_mat; beta <- beta[rownames(V)]
  tt <- delete.response(terms(f_fixed))
  a_mat <- t(sapply(group_levels, function(g) {
    md_g <- md; md_g$joint <- factor(g, levels = group_levels)
    colMeans(model.matrix(tt, md_g))[names(beta)]
  }))
  out <- lapply(names(ctr), function(nm) {
    d <- as.numeric(t(ctr[[nm]]) %*% a_mat)
    est <- sum(d * beta); se <- sqrt(as.numeric(t(d) %*% V %*% d))
    data.frame(contrast = nm, estimate = est, se = se,
               lci = est - 1.96 * se, uci = est + 1.96 * se,
               p = 2 * pnorm(-abs(est / se)), n = nrow(md))
  })
  do.call(rbind, out)
}

# climate PCs on the analytical sample (same convention as 01c)
climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")
cc <- complete.cases(ana[, climate_vars])
pc_fit <- prcomp(ana[cc, climate_vars], center = TRUE, scale. = TRUE)
ana[, paste0("climate_pc", 1:5)] <-
  predict(pc_fit, newdata = ana[, climate_vars])[, 1:5]
ana$ln_hh_consumption_total <- ifelse(ana$hh_consumption_total > 0,
                                      log(ana$hh_consumption_total), NA_real_)

# ---- 3. S18 --------------------------------------------------------------------
cat("\n=== S18: excluded vs included ===\n")
s4_w1 <- s4[s4$wave == 1 & s4$age > 45 & s4$age < 90, ]
excl <- s4_w1[rowSums(s4_w1[, diseases] == 1, na.rm = TRUE) > 0, ]
cat("Criterion-3 excluded at 2011 (biomarker-complete, 45-89):", nrow(excl), "\n")
excl <- merge(excl[, c("id", "age")], m11, by = "id", all.x = TRUE)
cat("Excluded with master covariates matched:", sum(!is.na(excl$bmi)), "\n")
excl$female <- as.integer(excl$sex_c == 0)   # master 性别: 0 = female

inc <- ana11
inc$female <- as.integer(inc$sex == 0)       # analysis dataset: 0 = female

# display-only BMI plausibility filter (12-60 kg/m2): descriptive tables
# exclude implausible values; propensity/outcome models use delivered values
inc$bmi_disp  <- ifelse(inc$bmi >= 12 & inc$bmi <= 60, inc$bmi, NA_real_)
excl$bmi_disp <- ifelse(excl$bmi >= 12 & excl$bmi <= 60, excl$bmi, NA_real_)

s18_specs <- list(
  list(label = "Age, years",                 ev = "age",             iv = "age",             bin = FALSE),
  list(label = "Female, n (%)",              ev = "female",          iv = "female",          bin = TRUE),
  list(label = "Body mass index, kg/m2 (plausible range 12-60)", ev = "bmi_disp", iv = "bmi_disp", bin = FALSE),
  list(label = "Current smoking, n (%)",     ev = "smoking_current", iv = "smoking_current", bin = TRUE),
  list(label = "Current alcohol use, n (%)", ev = "alcohol_current", iv = "alcohol_current", bin = TRUE),
  list(label = "Hypertension, n (%)",        ev = "hypertension",    iv = "hypertension",    bin = TRUE),
  list(label = "Dyslipidaemia, n (%)",       ev = "dyslipidemia",    iv = "dyslipidemia",    bin = TRUE),
  list(label = "EVI (1-km buffer, 4-year average)",     ev = "EVI_1KM_mean4y", iv = "EVI_1KM_mean4y", bin = FALSE),
  list(label = "O3 (1-km, 4-year average), ug/m3",      ev = "O3_1KM_mean4y",  iv = "O3_1KM_mean4y",  bin = FALSE))
s18 <- data.frame(Characteristic = character(), Excluded = character(),
                  Included = character(), Std_diff = numeric(),
                  N_excluded = integer(), N_included = integer())
for (sp in s18_specs) {
  r <- stat_row(excl[[sp$ev]], inc[[sp$iv]], sp$bin)
  s18 <- rbind(s18, data.frame(Characteristic = sp$label, Excluded = r$e,
                               Included = r$i, Std_diff = r$d,
                               N_excluded = r$n_e, N_included = r$n_i))
}
s18 <- rbind(s18, data.frame(
  Characteristic = "BAA (PhenoAge), years",
  Excluded = "not applicable (see note)",
  Included = sprintf("%.2f (%.2f)", mean(inc$BAA_PhenoAge, na.rm = TRUE),
                     sd(inc$BAA_PhenoAge, na.rm = TRUE)),
  Std_diff = NA_real_, N_excluded = NA_integer_,
  N_included = sum(!is.na(inc$BAA_PhenoAge))))
cat("\nS18 (N excluded =", nrow(excl), "; N included =", nrow(inc), "):\n")
print(s18, row.names = FALSE)
write.csv(s18, "output/tables/s18_excluded_vs_included_v4.csv",
          row.names = FALSE)

# ---- 4. attrition comparison (available-data basis) ----------------------------
cat("\n=== Attrition: retained vs 2011-only ===\n")
pool <- ref[ref$wave == 1 & ref$年龄 > 45 & ref$年龄 < 90,
            c("id", "年龄", "BAA_PhenoAge")]
names(pool) <- c("id", "age", "BAA_ref")
pool <- merge(pool, m11, by = "id")   # sex_c + covariates from the master
pool$ln_cons <- ifelse(pool$hh_consumption_total > 0,
                       log(pool$hh_consumption_total), NA_real_)
pool$retained <- pool$id %in% ana11$id
n_ret <- sum(pool$retained); n_att <- sum(!pool$retained)
cat("Eligible 2011 pool (post-QC, 45-89):", nrow(pool),
    " retained:", n_ret, " 2011-only:", n_att, "\n")

comp <- pool[pool$retained, ]
attrit <- pool[!pool$retained, ]
for (d in c("comp", "attrit")) {
  dd <- get(d)
  dd$female  <- as.integer(dd$sex_c == 0)
  dd$low_edu <- as.integer(dd$education == 1)
  dd$rural   <- as.integer(dd$urban_rural == 1)
  dd$bmi_disp <- ifelse(dd$bmi >= 12 & dd$bmi <= 60, dd$bmi, NA_real_)
  assign(d, dd)
}

attr_specs <- list(
  list(label = "Age, years",                 v = "age",             bin = FALSE),
  list(label = "Female, n (%)",              v = "female",          bin = TRUE),
  list(label = "Low education, n (%)",       v = "low_edu",         bin = TRUE),
  list(label = "Rural residence, n (%)",     v = "rural",           bin = TRUE),
  list(label = "Body mass index, kg/m2 (plausible range 12-60)", v = "bmi_disp", bin = FALSE),
  list(label = "Current smoking, n (%)",     v = "smoking_current", bin = TRUE),
  list(label = "Current alcohol use, n (%)", v = "alcohol_current", bin = TRUE),
  list(label = "Hypertension, n (%)",        v = "hypertension",    bin = TRUE),
  list(label = "Dyslipidaemia, n (%)",       v = "dyslipidemia",    bin = TRUE),
  list(label = "BAA (PhenoAge, recalibrated scale), years", v = "BAA_ref", bin = FALSE),
  list(label = "EVI (1-km buffer, 4-year average)",     v = "EVI_1KM_mean4y", bin = FALSE),
  list(label = "O3 (1-km, 4-year average), ug/m3",      v = "O3_1KM_mean4y",  bin = FALSE))
att <- data.frame(Characteristic = character(), Retained = character(),
                  Y2011_only = character(), Std_diff = numeric(),
                  N_retained = integer(), N_2011only = integer())
for (sp in attr_specs) {
  r <- stat_row(attrit[[sp$v]], comp[[sp$v]], sp$bin)
  att <- rbind(att, data.frame(Characteristic = sp$label, Retained = r$i,
                               Y2011_only = r$e, Std_diff = r$d,
                               N_retained = r$n_i, N_2011only = r$n_e))
}
cat("\nAttrition table:\n")
print(att, row.names = FALSE)
write.csv(att, "output/tables/attrition_baseline_comparison_v4.csv",
          row.names = FALSE)

# fill missing master mean4y exposures for retained participants with the
# analytical-dataset values (same participants, same construct; the master
# mean4y columns are incomplete for part of the analytical sample)
exp_fill <- ana11[, c("id", "EVI_1KM_mean4y", "O3_1KM_mean4y")]
names(exp_fill) <- c("id", "EVI_1KM_mean4y_ana", "O3_1KM_mean4y_ana")
pool <- merge(pool, exp_fill, by = "id", all.x = TRUE)
pool$EVI_1KM_mean4y <- ifelse(is.na(pool$EVI_1KM_mean4y),
                              pool$EVI_1KM_mean4y_ana, pool$EVI_1KM_mean4y)
pool$O3_1KM_mean4y <- ifelse(is.na(pool$O3_1KM_mean4y),
                             pool$O3_1KM_mean4y_ana, pool$O3_1KM_mean4y)
pool$EVI_1KM_mean4y_ana <- NULL; pool$O3_1KM_mean4y_ana <- NULL
# second tier: participants without a 4-year exposure history -> current-year
# values (propensity model and attriter joint-group classification only;
# the descriptive table above uses observed values only)
pool$EVI_1KM_mean4y <- ifelse(is.na(pool$EVI_1KM_mean4y),
                              pool$EVI_1KM, pool$EVI_1KM_mean4y)
pool$O3_1KM_mean4y <- ifelse(is.na(pool$O3_1KM_mean4y),
                             pool$O3_1KM, pool$O3_1KM_mean4y)

# ---- 5. IPAW joint models (primary-scale BAA) ----------------------------------
cat("\n=== IPAW joint models ===\n")
# parsimonious propensity model: excludes high-missingness covariates
# (pension insurance, household consumption, climate) so that nearly all
# retained participants receive a weight
ps_vars <- c("age", "sex_c", "education", "urban_rural", "marital_status",
             "smoking_current", "alcohol_current", "bmi", "hypertension",
             "dyslipidemia", "BAA_ref", "EVI_1KM_mean4y", "O3_1KM_mean4y")
psd <- pool[complete.cases(pool[, ps_vars]), ]
cat("Propensity-model frame:", nrow(psd), "of", nrow(pool),
    "(retained:", sum(psd$retained), ")\n")
ps_fit <- glm(retained ~ age + factor(sex_c) + factor(education) +
                factor(urban_rural) + factor(marital_status) +
                smoking_current + alcohol_current + bmi + hypertension +
                dyslipidemia + BAA_ref + EVI_1KM_mean4y + O3_1KM_mean4y,
              data = psd, family = binomial)
psd$ps <- predict(ps_fit, type = "response")
psd$ipaw <- 1 / psd$ps
qs <- quantile(psd$ipaw[psd$retained], c(0.01, 0.99))
psd$ipaw_w <- pmin(pmax(psd$ipaw, qs[1]), qs[2])
cat("IPAW (retained): mean", round(mean(psd$ipaw_w[psd$retained]), 3),
    " range", round(min(psd$ipaw_w[psd$retained]), 2), "-",
    round(max(psd$ipaw_w[psd$retained]), 2), "\n")

ana_ipw <- merge(ana, psd[psd$retained, c("id", "ipaw_w")], by = "id",
                 all.x = TRUE)
cat("Analytical rows with IPAW:", sum(!is.na(ana_ipw$ipaw_w)), "/",
    nrow(ana_ipw), "\n")

ipaw_rows <- list()
for (s in c("1KM", "500m")) {
  ana_ipw$w_comb <- ana_ipw[[paste0("ipw_joint_", s)]] * ana_ipw$ipaw_w
  res <- fit_joint(ana_ipw, paste0("joint_exp_", s), "w_comb")
  ipaw_rows[[s]] <- cbind(scale = s, res)
  cat("\nIPAW joint model", s, ":\n")
  print(ipaw_rows[[s]], row.names = FALSE)
}
ipaw_out <- do.call(rbind, ipaw_rows)
write.csv(ipaw_out, "output/tables/attrition_ipaw_joint_v4.csv",
          row.names = FALSE)

# ---- 6. pattern-mixture (recalibrated-BAA scale) -------------------------------
cat("\n=== Pattern-mixture (delta-adjustment, ref-BAA scale) ===\n")
# stayer frame: analytical rows with ref-scale BAA merged by (id, year)
ref_baa <- ref[, c("id", "wave", "BAA_PhenoAge")]
ref_baa$year <- ifelse(ref_baa$wave == 1, 2011, 2015)
stayer <- merge(ana, ref_baa[, c("id", "year", "BAA_PhenoAge")],
                by = c("id", "year"), suffixes = c("", "_refscale"))
stayer$BAA_ref <- stayer$BAA_PhenoAge_refscale
cat("Stayer rows with ref-scale BAA:", nrow(stayer), "\n")

# joint-group cutpoints from analytical baseline (mean4y exposures)
cut_evi <- median(ana11$EVI_1KM_mean4y, na.rm = TRUE)
cut_o3  <- median(ana11$O3_1KM_mean4y, na.rm = TRUE)
cut_evi5 <- median(ana11$EVI_500m_mean4y, na.rm = TRUE)
cut_o35  <- median(ana11$O3_500m_mean4y, na.rm = TRUE)
cat("Cutpoints 1KM: EVI", round(cut_evi, 4), " O3", round(cut_o3, 3),
    " 500m: EVI", round(cut_evi5, 4), " O3", round(cut_o35, 3), "\n")

# mean within-person BAA change among completers (MAR reference shift)
b11 <- stayer$BAA_ref[stayer$year == 2011]
b15 <- stayer$BAA_ref[stayer$year == 2015]
i11 <- stayer$id[stayer$year == 2011]; i15 <- stayer$id[stayer$year == 2015]
common <- intersect(i11, i15)
d_bar <- mean(b15[match(common, i15)] - b11[match(common, i11)], na.rm = TRUE)
cat("Mean within-person delta-BAA (ref scale) among completers:",
    round(d_bar, 3), "\n")

# climate PCs for pool (same loadings as the analytical sample)
cc2 <- complete.cases(pool[, climate_vars])
pool[, paste0("climate_pc", 1:5)] <- NA_real_
pool[cc2, paste0("climate_pc", 1:5)] <-
  predict(pc_fit, newdata = pool[cc2, climate_vars])[, 1:5]

mk_joint <- function(evi, o3, ce, co) {
  1L * (evi <= ce & o3 <= co) + 2L * (evi <= ce & o3 > co) +
  3L * (evi > ce & o3 <= co) + 4L * (evi > ce & o3 > co)
}
pool$joint_1KM <- mk_joint(pool$EVI_1KM_mean4y, pool$O3_1KM_mean4y,
                           cut_evi, cut_o3)
pool$joint_500m <- mk_joint(pool$EVI_500m_mean4y, pool$O3_500m_mean4y,
                            cut_evi5, cut_o35)

# attriters entering the augmented model need the full model covariate set
full_covs <- c("age", "sex_c", "education", "ln_cons", "retirement",
               "pension_insurance", "smoking_current", "alcohol_current",
               paste0("climate_pc", 1:5), "bmi", "hypertension",
               "dyslipidemia", "BAA_ref", "joint_1KM", "joint_500m")
attr_rows <- pool[!pool$retained &
                    complete.cases(pool[, full_covs]), ]
cat("Attriters with complete model covariates:", nrow(attr_rows), "of",
    sum(!pool$retained), "\n")

build_aug <- function(scale, delta) {
  jvar <- paste0("joint_exp_", scale)
  aj <- stayer[, c("BAA_ref", jvar, "id_num", "age", "sex", "education",
                   "ln_hh_consumption_total", "retirement",
                   "pension_insurance", "smoking_current", "alcohol_current",
                   paste0("climate_pc", 1:5), "bmi", "hypertension",
                   "dyslipidemia", "year")]
  names(aj)[1:2] <- c("BAA_use", "joint")
  a1 <- data.frame(
    BAA_use = attr_rows$BAA_ref,
    joint = attr_rows[[paste0("joint_", scale)]],
    id_num = attr_rows$id, age = attr_rows$age, sex = attr_rows$sex_c,
    education = attr_rows$education,
    ln_hh_consumption_total = attr_rows$ln_cons,
    retirement = attr_rows$retirement,
    pension_insurance = attr_rows$pension_insurance,
    smoking_current = attr_rows$smoking_current,
    alcohol_current = attr_rows$alcohol_current,
    attr_rows[, paste0("climate_pc", 1:5)],
    bmi = attr_rows$bmi, hypertension = attr_rows$hypertension,
    dyslipidemia = attr_rows$dyslipidemia, year = 2011)
  a2 <- a1
  a2$BAA_use <- attr_rows$BAA_ref + d_bar + delta   # delta-adjusted 2015 BAA
  a2$age <- attr_rows$age + 4
  a2$year <- 2015
  rbind(aj, a1, a2)
}

pm_rows <- list()
for (s in c("1KM", "500m")) {
  for (dl in 0:3) {
    aug <- build_aug(s, dl)
    aug$w <- 1
    res <- fit_joint(aug, "joint", "w", outcome = "BAA_use")
    did <- res[res$contrast == "Difference-in-differences interaction", ]
    pm_rows[[length(pm_rows) + 1]] <- data.frame(
      scale = s, delta = dl, did_estimate = did$estimate,
      did_lci = did$lci, did_uci = did$uci, did_p = did$p, n = did$n)
    cat(sprintf("scale=%s delta=%d  DiD=%.3f [%.3f, %.3f] p=%.3f n=%d\n",
                s, dl, did$estimate, did$lci, did$uci, did$p,
                as.integer(did$n)))
  }
}
pm_out <- do.call(rbind, pm_rows)
write.csv(pm_out, "output/tables/attrition_pattern_mixture_v4.csv",
          row.names = FALSE)

cat("\nAll attrition outputs written.\n")

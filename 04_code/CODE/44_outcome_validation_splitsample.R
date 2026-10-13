# Full9 to reduced7 excludes immune_balance and rdw_synthetic only.
# The albumin proxy remains in both specifications; counts of rows are observations.
# =============================================================================
# 44_outcome_validation_splitsample.R  (v5 revision, task T1)
# Auxiliary split-sample calibration of a study-specific biomarker score.
# The original calibration of the stored primary BAA is unavailable; this
# procedure is not verified as a reconstruction or validation of that outcome.
# This script fits an elastic-net Cox model (alpha=0.5, 5-fold CV):
#   - 2/3 train / 1/3 test split clustered by participant id
#   - recalibrate on train, evaluate discrimination on held-out test
#   - alpha sensitivity {0.1, 0.5, 0.9}; reduced biomarker set (drops the two
#     constructed substitutes: immune_balance, rdw_synthetic)
#   - event-anchored HR per SD BAA (age-adjusted) on the test set
# Raw CHARLS biomarker file is restricted: the script SKIPs gracefully if absent.
# =============================================================================

suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr); library(survival); library(glmnet)
})

set.seed(20240101)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/44_outcome_validation_splitsample.log", split = TRUE)
on.exit(sink(), add = TRUE)

RAW <- "DATA/restricted/biomarker_selected.dta"
if (!file.exists(RAW)) {
  cat("SKIP: restricted raw biomarker file not found at", RAW, "\n")
  quit(save = "no", status = 0)
}

raw <- as.data.frame(read_dta(RAW))
cat("Raw rows:", nrow(raw), "\n")

# ---- implement the documented auxiliary transforms and filters ----
d <- raw %>%
  mutate(
    exam_year = case_when(wave == 1 ~ 2011, wave == 3 ~ 2015, TRUE ~ 2015),
    albumin = bl_cysc,
    creatinine = bl_crea * 88.4 * 0.85,
    glucose = bl_glu / 18,
    log_crp = log(bl_crp / 10 + 1e-5),
    wbc = bl_wbc, mcv = bl_mcv, hgb = bl_hgb,
    rdw_synthetic = (mcv / hgb) * 100,
    immune_balance = wbc / (bl_crp + 1),
    death_status = ifelse(is.na(radyear), 0, 1),
    death_time = case_when(
      death_status == 1 & radyear > exam_year ~ radyear - exam_year,
      death_status == 1 & radyear <= exam_year ~ NA_real_,
      death_status == 0 ~ 2020 - exam_year,
      TRUE ~ NA_real_)) %>%
  filter(!is.na(age), bl_fasting == 1,
         between(albumin, 0.5, 4.0), bl_crp < 10,
         between(creatinine, 30, 300),
         !is.na(death_time), death_time > 0,
         !is.na(exam_year))

cat("Filtered calibration frame:", nrow(d), "rows | deaths:",
    sum(d$death_status == 1), "\n")

full_set <- c("albumin", "creatinine", "glucose", "log_crp", "immune_balance",
              "mcv", "rdw_synthetic", "wbc", "age")
reduced_set <- setdiff(full_set, c("immune_balance", "rdw_synthetic"))

# ---- clustered 2/3 - 1/3 split by participant id ----
ids <- unique(d$id)
train_ids <- sample(ids, size = floor(length(ids) * 2 / 3))
d$split <- ifelse(d$id %in% train_ids, "train", "test")
cat("Split: train", sum(d$split == "train"), "rows | test",
    sum(d$split == "test"), "rows | test deaths:",
    sum(d$death_status[d$split == "test"] == 1), "\n")

run_variant <- function(preds, alpha, label) {
  tr <- d[d$split == "train", c(preds, "death_time", "death_status")] %>% drop_na()
  te <- d[d$split == "test",  c(preds, "death_time", "death_status")] %>% drop_na()
  Xtr <- as.matrix(tr[, preds]); Xte <- as.matrix(te[, preds])
  cv <- cv.glmnet(Xtr, Surv(tr$death_time, tr$death_status),
                  family = "cox", alpha = alpha, nfolds = 5,
                  maxit = 1000, thresh = 1e-6)
  lp_tr <- as.numeric(predict(cv, Xtr, s = "lambda.min"))
  lp_te <- as.numeric(predict(cv, Xte, s = "lambda.min"))
  baa_te <- (lp_te - mean(lp_tr)) / sd(lp_tr) * 15   # BAA scale, train moments
  # age is among predictors; use raw age column for the age-only reference
  age_te <- te$age
  s_te <- Surv(te$death_time, te$death_status)
  c_baa  <- summary(coxph(s_te ~ baa_te))$concordance[1]
  c_age  <- summary(coxph(s_te ~ age_te))$concordance[1]
  c_both <- summary(coxph(s_te ~ baa_te + age_te))$concordance[1]
  hr <- summary(coxph(s_te ~ scale(baa_te) + age_te))$coefficients
  data.frame(
    variant = label, alpha = alpha, n_predictors = length(preds),
    n_train = nrow(tr), n_test = nrow(te), deaths_test = sum(te$death_status),
    c_baa_alone = c_baa, c_age_alone = c_age, c_baa_plus_age = c_both,
    delta_c_over_age = c_both - c_age,
    hr_per_sd_baa = hr[1, "coef"] |> exp(),
    hr_lci = exp(hr[1, "coef"] - 1.96 * hr[1, "se(coef)"]),
    hr_uci = exp(hr[1, "coef"] + 1.96 * hr[1, "se(coef)"]),
    hr_p = hr[1, "Pr(>|z|)"])
}

res <- rbind(
  run_variant(full_set,    0.5, "Full biomarker set, alpha=0.5"),
  run_variant(full_set,    0.1, "Full biomarker set, alpha=0.1"),
  run_variant(full_set,    0.9, "Full biomarker set, alpha=0.9"),
  run_variant(reduced_set, 0.5, "Reduced set (immune_balance and rdw_synthetic removed), alpha=0.5"))

print(res, digits = 4)
write.csv(res, "output/tables/outcome_validation_splitsample_v5.csv",
          row.names = FALSE)
cat("\nDone. output/tables/outcome_validation_splitsample_v5.csv\n")

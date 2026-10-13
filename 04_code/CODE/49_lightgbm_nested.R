# Validation boundary: the supplied single imputation and initial climate PCA
# precede participant splitting. This is not fully leakage-free validation.
# The filename is retained for compatibility; this is not nested cross-validation.
# =============================================================================
# 49_lightgbm_nested.R  (v5 revision, task T5)
# Participant-level split re-evaluation of the LightGBM BAA-prediction model:
#   - participant-level 60/20/20 train/validation/test split (repeated
#     measures of the same participant never cross splits)
#   - early stopping on the VALIDATION set (test set untouched during tuning)
#   - honest held-out performance (R2/RMSE/MAE) on the test set
#   - permutation importance on the test set (10 repeats)
#   - SHAP values on the test set (exact TreeSHAP, type = "contrib"),
#     aggregated into feature groups
#   - counterfactual scenario (EVI >= 0.40, O3 <= 60 ug/m3) on the test set
#     with 1,000-rep participant-level bootstrap CI
# Output: output/tables/lightgbm_nested_*.csv
# =============================================================================

suppressPackageStartupMessages({
  library(haven); library(dplyr); library(lightgbm)
})

set.seed(20240101)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/49_lightgbm_nested.log", split = TRUE)
on.exit(sink(), add = TRUE)

# ---- data (same preprocessing as 34_ml_counterfactual_shap.R) ----
data <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
data <- data[data$age > 45 & data$age < 90, ]
data$ln_hh_consumption_total <- ifelse(data$hh_consumption_total <= 0, NA,
                                       log(data$hh_consumption_total))
cv <- c("TMP", "PER", "SNOW", "DEM", "lat")
pc <- prcomp(data[complete.cases(data[, cv]), cv], center = TRUE, scale. = TRUE)
data[, paste0("climate_pc", 1:5)] <- predict(pc, newdata = data[, cv])[, 1:5]

vars <- c("BAA_PhenoAge", "EVI_1KM_mean4y", "O3_1KM_mean4y",
          "age", "sex", "ln_hh_consumption_total", "retirement",
          "pension_insurance", "smoking_current", "alcohol_current",
          paste0("climate_pc", 1:5), "bmi", "hypertension", "dyslipidemia",
          "year", "id_num")
d <- data[, vars]
d <- d[complete.cases(d), ]
cat("Analysis frame:", nrow(d), "obs /", length(unique(d$id_num)), "participants\n")

fac <- c("sex", "retirement", "pension_insurance", "smoking_current",
         "alcohol_current", "hypertension", "dyslipidemia", "year")
for (v in fac) d[[v]] <- as.integer(as.factor(d[[v]])) - 1

# ---- participant-level split 60/20/20 ----
ids <- unique(d$id_num)
n_ids <- length(ids)
train_ids <- sample(ids, floor(0.6 * n_ids))
rest <- setdiff(ids, train_ids)
val_ids <- sample(rest, floor(0.5 * length(rest)))
test_ids <- setdiff(rest, val_ids)
tr <- d[d$id_num %in% train_ids, ]
va <- d[d$id_num %in% val_ids, ]
te <- d[d$id_num %in% test_ids, ]
cat("Split (obs): train", nrow(tr), "/ val", nrow(va), "/ test", nrow(te),
    "| participants:", length(train_ids), "/", length(val_ids), "/",
    length(test_ids), "\n")

X <- function(x) as.matrix(x[, setdiff(names(x), c("BAA_PhenoAge", "id_num"))])
dtrain <- lgb.Dataset(X(tr), label = tr$BAA_PhenoAge)
dval   <- lgb.Dataset.create.valid(dtrain, X(va), label = va$BAA_PhenoAge)

params <- list(objective = "regression", metric = "l2", learning_rate = 0.05,
               num_leaves = 31, max_depth = -1, min_data_in_leaf = 20,
               feature_fraction = 0.8, bagging_fraction = 0.8, bagging_freq = 5,
               lambda_l1 = 0.1, lambda_l2 = 0.1)
fit <- lgb.train(params, dtrain, valids = list(validation = dval),
                 nrounds = 1000, early_stopping_rounds = 50, verbose = -1)
cat("Best iteration (validation early stopping):", fit$best_iter, "\n")

pred_te <- predict(fit, X(te))
rmse <- sqrt(mean((te$BAA_PhenoAge - pred_te)^2))
mae  <- mean(abs(te$BAA_PhenoAge - pred_te))
r2   <- 1 - sum((te$BAA_PhenoAge - pred_te)^2) /
            sum((te$BAA_PhenoAge - mean(te$BAA_PhenoAge))^2)
perf <- data.frame(model = "LightGBM nested (participant-level split)",
                   n_train = nrow(tr), n_val = nrow(va), n_test = nrow(te),
                   best_iter = fit$best_iter, r2_test = r2, rmse_test = rmse,
                   mae_test = mae)
print(perf, digits = 4)
write.csv(perf, "output/tables/lightgbm_nested_performance_v5.csv",
          row.names = FALSE)

# ---- permutation importance on test set (10 repeats) ----
feat <- colnames(X(te))
perm <- lapply(feat, function(f) {
  drops <- replicate(10, {
    Xp <- X(te); Xp[, f] <- sample(Xp[, f])
    pp <- predict(fit, Xp)
    r2p <- 1 - sum((te$BAA_PhenoAge - pp)^2) /
               sum((te$BAA_PhenoAge - mean(te$BAA_PhenoAge))^2)
    r2 - r2p
  })
  data.frame(feature = f, r2_drop_mean = mean(drops), r2_drop_sd = sd(drops))
})
perm <- do.call(rbind, perm)
perm <- perm[order(-perm$r2_drop_mean), ]
cat("\nPermutation importance (test-set R2 drop):\n"); print(perm, digits = 4)
write.csv(perm, "output/tables/lightgbm_nested_importance_v5.csv",
          row.names = FALSE)

# ---- SHAP on test set (exact TreeSHAP via lightgbm predcontrib) + groups ----
cat("\nComputing exact TreeSHAP contributions on test set...\n")
sv_mat <- predict(fit, X(te), type = "contrib")
sv <- sv_mat[, -ncol(sv_mat), drop = FALSE]   # drop bias column
colnames(sv) <- colnames(X(te))
mas <- sort(colMeans(abs(sv)), decreasing = TRUE)
groups <- list(
  "Exposures (EVI + O3)"      = c("EVI_1KM_mean4y", "O3_1KM_mean4y"),
  "Climate (PC1-PC5)"         = paste0("climate_pc", 1:5),
  "Demographics (age, sex)"   = c("age", "sex"),
  "Socioeconomic"             = c("ln_hh_consumption_total", "retirement",
                                  "pension_insurance"),
  "Lifestyle (smoking, alcohol)" = c("smoking_current", "alcohol_current"),
  "Health status"             = c("bmi", "hypertension", "dyslipidemia"),
  "Survey year"               = "year")
gshap <- do.call(rbind, lapply(names(groups), function(g) {
  fs <- intersect(groups[[g]], names(mas))
  data.frame(group = g, mean_abs_shap = sum(mas[fs]),
             share_pct = 100 * sum(mas[fs]) / sum(mas),
             features = paste(fs, collapse = ", "))
}))
gshap <- gshap[order(-gshap$mean_abs_shap), ]
feat_shap <- data.frame(feature = names(mas), mean_abs_shap = as.numeric(mas),
                        share_pct = 100 * as.numeric(mas) / sum(mas))
cat("\nGrouped SHAP:\n"); print(gshap, digits = 4)
write.csv(gshap, "output/tables/lightgbm_nested_shap_groups_v5.csv",
          row.names = FALSE)
write.csv(feat_shap, "output/tables/lightgbm_nested_shap_features_v5.csv",
          row.names = FALSE)

# ---- counterfactual scenario on test set + participant bootstrap CI ----
cf <- te
cf$EVI_1KM_mean4y <- pmax(cf$EVI_1KM_mean4y, 0.40)
cf$O3_1KM_mean4y  <- pmin(cf$O3_1KM_mean4y, 60)
pred_cf <- predict(fit, X(cf))
red <- data.frame(id = te$id_num, reduction = pred_te - pred_cf)
red <- red[order(red$id), ]
mean_red <- mean(red$reduction)

B <- 1000
tid <- unique(red$id)
boot <- replicate(B, {
  b_ids <- sample(tid, replace = TRUE)
  idx <- unlist(lapply(b_ids, function(i) which(red$id == i)))
  mean(red$reduction[idx])
})
ci <- quantile(boot, c(0.025, 0.975))
cf_out <- data.frame(
  scenario = "EVI >= 0.40 and O3 <= 60 ug/m3 (test set)",
  n_test = nrow(te), mean_reduction = mean_red,
  ci_lower = ci[1], ci_upper = ci[2],
  pct_improved = 100 * mean(red$reduction > 0))
cat("\nCounterfactual (test set):\n"); print(cf_out, digits = 4)
write.csv(cf_out, "output/tables/lightgbm_nested_counterfactual_v5.csv",
          row.names = FALSE)

cat("\nDone. output/tables/lightgbm_nested_*_v5.csv\n")

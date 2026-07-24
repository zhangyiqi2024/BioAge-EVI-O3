# =============================================================================
# v4 predictive scenario contrasts for BAA
# Purpose: train a predictive LightGBM model and report descriptive model-based
#          scenario contrasts. These are not causal counterfactual estimates.
# =============================================================================

library(haven)
library(dplyr)
library(caret)
library(lightgbm)
library(boot)

set.seed(20240101)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

sink("logs/36_predictive_scenarios.log", split = TRUE)
on.exit(sink(), add = TRUE)

raw <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
raw$ln_hh_consumption_total <- ifelse(raw$hh_consumption_total > 0,
                                      log(raw$hh_consumption_total), NA_real_)

climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")
pc_fit <- prcomp(raw[complete.cases(raw[, climate_vars]), climate_vars], center = TRUE, scale. = TRUE)
pc_scores <- predict(pc_fit, newdata = raw[, climate_vars])
raw <- cbind(raw, as.data.frame(pc_scores[, 1:5, drop = FALSE]))
names(raw)[(ncol(raw) - 4):ncol(raw)] <- paste0("climate_pc", 1:5)

features <- c("EVI_1KM_mean4y", "O3_1KM_mean4y", "age", "sex",
              "ln_hh_consumption_total", "retirement", "pension_insurance",
              "smoking_current", "alcohol_current", paste0("climate_pc", 1:5),
              "bmi", "hypertension", "dyslipidemia", "year")
keep <- c("BAA_PhenoAge", features)
dat <- raw[, keep]
dat <- dat[complete.cases(dat), ]
# Apply the cohort age-inclusion criterion (46-89 years), consistent with
# all other analyses; previously omitted, which admitted 88 out-of-range obs.
dat <- dat[dat$age > 45 & dat$age < 90, ]

factor_vars <- c("sex", "retirement", "pension_insurance", "smoking_current",
                 "alcohol_current", "hypertension", "dyslipidemia", "year")
for (v in factor_vars) dat[[v]] <- as.integer(as.factor(dat[[v]])) - 1L

set.seed(20240101)
train_id <- createDataPartition(dat$BAA_PhenoAge, p = 0.80, list = FALSE)
train <- dat[train_id, ]
test <- dat[-train_id, ]

x_train <- as.matrix(dplyr::select(train, -BAA_PhenoAge))
x_test <- as.matrix(dplyr::select(test, -BAA_PhenoAge))
dtrain <- lgb.Dataset(x_train, label = train$BAA_PhenoAge)
dtest <- lgb.Dataset.create.valid(dtrain, x_test, label = test$BAA_PhenoAge)

params <- list(objective = "regression", metric = "l2", learning_rate = 0.05,
               num_leaves = 31, min_data_in_leaf = 20, feature_fraction = 0.8,
               bagging_fraction = 0.8, bagging_freq = 5, lambda_l1 = 0.1,
               lambda_l2 = 0.1)
model <- lgb.train(params = params, data = dtrain, valids = list(test = dtest),
                   nrounds = 1000, early_stopping_rounds = 50, verbose = -1)

predict_scenario <- function(data, set_evi = FALSE, set_o3 = FALSE,
                             target_evi = 0.40, target_o3 = 60) {
  cf <- data
  eligible_evi <- cf$EVI_1KM_mean4y < target_evi
  eligible_o3 <- cf$O3_1KM_mean4y > target_o3
  if (set_evi) cf$EVI_1KM_mean4y[eligible_evi] <- target_evi
  if (set_o3) cf$O3_1KM_mean4y[eligible_o3] <- target_o3
  pred_current <- predict(model, as.matrix(dplyr::select(data, -BAA_PhenoAge)))
  pred_scenario <- predict(model, as.matrix(dplyr::select(cf, -BAA_PhenoAge)))
  data.frame(
    current_prediction = pred_current,
    scenario_prediction = pred_scenario,
    predicted_difference = pred_scenario - pred_current,
    eligible_evi = eligible_evi,
    eligible_o3 = eligible_o3
  )
}

scenarios <- list(
  "EVI-only: values below 0.40 set to 0.40" = predict_scenario(dat, set_evi = TRUE),
  "O3-only: values above 60 set to 60" = predict_scenario(dat, set_o3 = TRUE),
  "Joint: EVI and O3 targets applied" = predict_scenario(dat, set_evi = TRUE, set_o3 = TRUE)
)

summary_rows <- lapply(names(scenarios), function(name) {
  z <- scenarios[[name]]
  n <- nrow(z)
  mean_diff <- mean(z$predicted_difference)
  se_diff <- sd(z$predicted_difference) / sqrt(n)
  data.frame(
    scenario = name,
    N = n,
    eligible_evi_n = sum(z$eligible_evi),
    eligible_o3_n = sum(z$eligible_o3),
    eligible_both_n = sum(z$eligible_evi & z$eligible_o3),
    mean_predicted_BAA_difference = mean_diff,
    sd_predicted_BAA_difference = sd(z$predicted_difference),
    median_predicted_BAA_difference = median(z$predicted_difference),
    ci_lower_naive = mean_diff - 1.96 * se_diff,
    ci_upper_naive = mean_diff + 1.96 * se_diff
  )
})
summary_table <- bind_rows(summary_rows)
write.csv(summary_table, "output/tables/predictive_scenario_contrasts_v4.csv", row.names = FALSE)

performance <- data.frame(
  N_train = nrow(train), N_test = nrow(test), best_iteration = model$best_iter,
  RMSE_test = sqrt(mean((predict(model, x_test) - test$BAA_PhenoAge)^2)),
  R2_test = 1 - sum((predict(model, x_test) - test$BAA_PhenoAge)^2) /
    sum((test$BAA_PhenoAge - mean(test$BAA_PhenoAge))^2)
)
write.csv(performance, "output/tables/predictive_model_performance_v4.csv", row.names = FALSE)

cat("Predictive scenario analysis completed.\n")
print(summary_table)
print(performance)

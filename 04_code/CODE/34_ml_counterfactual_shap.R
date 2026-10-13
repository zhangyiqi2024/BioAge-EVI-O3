############################################################
#  Machine Learning Counterfactual Prediction + Climate PC + LightGBM + SHAP
#  Revised version & supplementary figures (incl. Table X)
############################################################

#===================== 0. Working directory & packages =====================#

# setwd("DATA/analysis_ready/")  # Working directory should be set to project root

library(haven)        # read .dta
library(glmnet)       # LASSO / ENet
library(gbm)          # GBM
library(pls)          # PLS
library(Cubist)       # Cubist
library(mgcv)         # GAM
library(lightgbm)     # LightGBM
library(dplyr)
library(ggplot2)
library(caret)
library(patchwork)
library(viridis)
library(fastshap)
library(MLmetrics)
library(scales)
library(grid)
library(tidyr)        # pivot_longer etc.

set.seed(20240101)

#===================== 1. Read data & revised PCA =====================#

cat("Reading data...\n")
data <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
cat("Data read complete, sample size:", nrow(data), "variables:", ncol(data), "\n")

# Apply the cohort age-inclusion criterion (46-89 years), consistent with
# all other analyses (see 36_predictive_scenarios.R).
data <- data[data$age > 45 & data$age < 90, ]
cat("After age-inclusion filter (46-89 years):", nrow(data), "observations\n")

# list.files()  # commented out — use project-relative data paths

# Log-transformed household consumption
if ("hh_consumption_total" %in% names(data)) {
  data$ln_hh_consumption_total <- ifelse(data$hh_consumption_total <= 0, NA, log(data$hh_consumption_total))
  cat("Created log household consumption variable: ln_hh_consumption_total\n")
}

# -------- Fix 1: PCA scores aligned to each individual -------- #
cat("Running PCA on climate variables...\n")
climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")
existing_climate_vars <- climate_vars[climate_vars %in% names(data)]

if (length(existing_climate_vars) > 0) {
  # Fit PCA on complete cases
  climate_complete <- data[complete.cases(data[, existing_climate_vars]),
                           existing_climate_vars, drop = FALSE]

  pca_result <- prcomp(climate_complete, scale. = TRUE, center = TRUE)

  # Project all samples so row count matches data exactly
  climate_all <- data[, existing_climate_vars, drop = FALSE]
  climate_pc_scores <- predict(pca_result, newdata = climate_all)
  climate_pc_scores <- as.data.frame(climate_pc_scores)
  colnames(climate_pc_scores) <- paste0("climate_pc", seq_len(ncol(climate_pc_scores)))

  data <- cbind(data, climate_pc_scores)

  pca_summary <- summary(pca_result)
  cat("Climate PCA complete\n")
  cat("Variance explained by first 5 PCs:\n")
  print(pca_summary$importance[2, 1:min(5, ncol(pca_result$x))])
} else {
  cat("Warning: climate variables not found, skipping PCA\n")
}

#===================== 2. Preprocessing function (mixed types) =====================#

preprocess_ml_data_mixed <- function(data) {
  # Variable selection (excluding hobbies, including climate PCs)
  selected_vars <- c(
    "BAA_PhenoAge", "EVI_1KM_mean4y", "O3_1KM_mean4y",
    "age", "sex", "ln_hh_consumption_total", "retirement", "pension_insurance",
    "smoking_current", "alcohol_current",
    "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4", "climate_pc5",
    "bmi", "hypertension", "dyslipidemia", "year"
  )

  existing_vars <- selected_vars[selected_vars %in% names(data)]
  cat("Variables used:", paste(existing_vars, collapse = ", "), "\n")

  # -------- Fix 2: base subset instead of select(all_of(...)) -------- #
  analysis_data <- data[, existing_vars, drop = FALSE]
  analysis_data <- na.omit(analysis_data)

  cat("Analysis sample size:", nrow(analysis_data), "\n")

  numeric_data <- analysis_data

  factor_vars <- c("sex", "retirement", "pension_insurance", "smoking_current",
                   "alcohol_current", "hypertension", "dyslipidemia", "year")
  factor_vars <- factor_vars[factor_vars %in% names(numeric_data)]

  for (var in factor_vars) {
    numeric_data[[var]] <- as.numeric(as.factor(numeric_data[[var]]))
  }

  factor_data <- analysis_data
  for (var in factor_vars) {
    factor_data[[var]] <- as.factor(factor_data[[var]])
  }

  lightgbm_data <- analysis_data
  for (var in factor_vars) {
    lightgbm_data[[var]] <- as.integer(as.factor(lightgbm_data[[var]])) - 1
  }

  list(
    numeric_data = numeric_data,
    factor_data = factor_data,
    lightgbm_data = lightgbm_data,
    factor_vars = factor_vars
  )
}

#===================== 3. Comprehensive performance metrics =====================#

calculate_comprehensive_metrics <- function(y_true, y_pred) {
  rmse <- sqrt(mean((y_true - y_pred)^2))
  mae  <- mean(abs(y_true - y_pred))
  r2   <- 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2)
  
  mape <- mean(abs((y_true - y_pred) / y_true)) * 100
  rae  <- sum(abs(y_true - y_pred)) / sum(abs(y_true - mean(y_true)))
  rrse <- sqrt(sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2))
  
  n <- length(y_true); p <- 1
  adj_r2 <- 1 - (1 - r2) * (n - 1) / (n - p - 1)
  explained_variance <- 1 - var(y_true - y_pred) / var(y_true)
  
  c(
    RMSE = rmse,
    MAE  = mae,
    R2   = r2,
    Adj_R2 = adj_r2,
    Explained_Variance = explained_variance,
    MAPE = mape,
    RAE  = rae,
    RRSE = rrse
  )
}

#===================== 4. Model training function =====================#

train_ml_models_mixed <- function(processed_data) {
  numeric_data   <- processed_data$numeric_data
  factor_data    <- processed_data$factor_data
  lightgbm_data  <- processed_data$lightgbm_data
  
  target <- numeric_data$BAA_PhenoAge
  train_indices <- createDataPartition(target, p = 0.8, list = FALSE)
  
  train_numeric  <- numeric_data[train_indices, ]
  test_numeric   <- numeric_data[-train_indices, ]
  train_factor   <- factor_data[train_indices, ]
  test_factor    <- factor_data[-train_indices, ]
  train_lightgbm <- lightgbm_data[train_indices, ]
  test_lightgbm  <- lightgbm_data[-train_indices, ]
  
  cat("Training set size:", nrow(train_numeric), "\n")
  cat("Test set size:", nrow(test_numeric),  "\n")
  
  x_train_numeric <- as.matrix(train_numeric %>% dplyr::select(-BAA_PhenoAge))
  x_test_numeric  <- as.matrix(test_numeric  %>% dplyr::select(-BAA_PhenoAge))
  y_train <- train_numeric$BAA_PhenoAge
  y_test  <- test_numeric$BAA_PhenoAge
  
  models_list      <- list()
  performance_list <- list()
  
  # 1. LASSO ------------------------------------------------
  cat("1. Training LASSO model...\n")
  tryCatch({
    lasso_model <- cv.glmnet(x_train_numeric, y_train,
                             alpha = 1, family = "gaussian")
    lasso_pred <- predict(lasso_model, newx = x_test_numeric, s = "lambda.min")
    lasso_metrics <- calculate_comprehensive_metrics(y_test, as.numeric(lasso_pred))
    models_list$lasso <- lasso_model
    performance_list$LASSO <- lasso_metrics
  }, error = function(e) cat("LASSO training failed:", e$message, "\n"))
  
  # 2. Elastic Net ------------------------------------------
  cat("2. Training Elastic Net model...\n")
  tryCatch({
    enet_model <- cv.glmnet(x_train_numeric, y_train,
                            alpha = 0.5, family = "gaussian")
    enet_pred <- predict(enet_model, newx = x_test_numeric, s = "lambda.min")
    enet_metrics <- calculate_comprehensive_metrics(y_test, as.numeric(enet_pred))
    models_list$enet <- enet_model
    performance_list$ElasticNet <- enet_metrics
  }, error = function(e) cat("Elastic Net training failed:", e$message, "\n"))
  
  # 3. GBM ---------------------------------------------------
  cat("3. Training GBM model...\n")
  tryCatch({
    gbm_model <- gbm(
      BAA_PhenoAge ~ .,
      data = train_factor,
      distribution = "gaussian",
      n.trees = 1000,
      interaction.depth = 4,
      shrinkage = 0.01,
      cv.folds = 5,
      n.minobsinnode = 10,
      verbose = FALSE
    )
    best_trees <- gbm.perf(gbm_model, method = "cv", plot.it = FALSE)
    gbm_pred <- predict(gbm_model, test_factor, n.trees = best_trees)
    gbm_metrics <- calculate_comprehensive_metrics(y_test, as.numeric(gbm_pred))
    models_list$gbm <- gbm_model
    performance_list$GBM <- gbm_metrics
    cat("GBM training succeeded, optimal trees:", best_trees, "\n")
  }, error = function(e) cat("GBM training failed:", e$message, "\n"))
  
  # 4. LightGBM ----------------------------------------------
  cat("4. Training LightGBM model...\n")
  tryCatch({
    lgb_train_data <- as.matrix(train_lightgbm %>% dplyr::select(-BAA_PhenoAge))
    lgb_test_data  <- as.matrix(test_lightgbm  %>% dplyr::select(-BAA_PhenoAge))
    dtrain <- lgb.Dataset(data = lgb_train_data, label = y_train)
    dtest  <- lgb.Dataset.create.valid(dtrain, data = lgb_test_data, label = y_test)
    
    params <- list(
      objective = "regression",
      metric = "l2",
      learning_rate = 0.05,
      num_leaves = 31,
      max_depth = -1,
      min_data_in_leaf = 20,
      feature_fraction = 0.8,
      bagging_fraction = 0.8,
      bagging_freq = 5,
      lambda_l1 = 0.1,
      lambda_l2 = 0.1
    )
    
    lgb_model <- lgb.train(
      params = params,
      data = dtrain,
      valids = list(test = dtest),
      nrounds = 1000,
      early_stopping_rounds = 50,
      verbose = -1
    )
    
    lgb_pred <- predict(lgb_model, lgb_test_data)
    lgb_metrics <- calculate_comprehensive_metrics(y_test, as.numeric(lgb_pred))
    models_list$lightgbm <- lgb_model
    performance_list$LightGBM <- lgb_metrics
    
    cat("LightGBM optimal iterations:", lgb_model$best_iter, "\n")
  }, error = function(e) cat("LightGBM training failed:", e$message, "\n"))
  
  # 5. Cubist ------------------------------------------------
  cat("5. Training Cubist model...\n")
  tryCatch({
    cubist_model <- cubist(x = x_train_numeric, y = y_train, committees = 10)
    cubist_pred  <- predict(cubist_model, x_test_numeric)
    cubist_metrics <- calculate_comprehensive_metrics(y_test, as.numeric(cubist_pred))
    models_list$cubist <- cubist_model
    performance_list$Cubist <- cubist_metrics
  }, error = function(e) cat("Cubist training failed:", e$message, "\n"))
  
  # 6. PLS ---------------------------------------------------
  cat("6. Training PLS model...\n")
  tryCatch({
    pls_model <- plsr(BAA_PhenoAge ~ ., data = train_numeric,
                      validation = "CV", ncomp = 10)
    best_comp <- which.min(RMSEP(pls_model)$val[1,,]) - 1
    pls_pred <- predict(pls_model, test_numeric, ncomp = best_comp)
    pls_metrics <- calculate_comprehensive_metrics(y_test, as.numeric(pls_pred))
    models_list$pls <- pls_model
    performance_list$PLS <- pls_metrics
  }, error = function(e) cat("PLS training failed:", e$message, "\n"))
  
  # 7. GAM ---------------------------------------------------
  cat("7. Training GAM model...\n")
  tryCatch({
    continuous_vars <- c(
      "EVI_1KM_mean4y", "O3_1KM_mean4y", "age", "ln_hh_consumption_total",
      "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4", "climate_pc5",
      "bmi"
    )
    smooth_terms <- paste0("s(", continuous_vars, ")", collapse = " + ")
    linear_terms <- paste(
      setdiff(names(train_numeric), c("BAA_PhenoAge", continuous_vars)),
      collapse = " + "
    )
    gam_formula <- as.formula(paste("BAA_PhenoAge ~", paste(smooth_terms, linear_terms, sep = " + ")))
    gam_model <- gam(gam_formula, data = train_numeric, method = "REML")
    gam_pred  <- predict(gam_model, newdata = test_numeric)
    gam_metrics <- calculate_comprehensive_metrics(y_test, as.numeric(gam_pred))
    models_list$gam <- gam_model
    performance_list$GAM <- gam_metrics
  }, error = function(e) cat("GAM training failed:", e$message, "\n"))
  
  # Summarize performance
  performance_df <- do.call(rbind, performance_list)
  performance_df <- as.data.frame(performance_df)
  performance_df$Model <- rownames(performance_df)
  rownames(performance_df) <- NULL
  
  cat("\nModel performance comparison:\n")
  print(performance_df)
  
  best_model_idx  <- which.max(performance_df$R2)
  best_model_name <- performance_df$Model[best_model_idx]
  best_model_r2   <- performance_df$R2[best_model_idx]
  cat("\nBest model:", best_model_name, "(R² =", round(best_model_r2, 4), ")\n")
  
  list(
    models = models_list,
    best_model_name = best_model_name,
    performance = performance_df,
    train_numeric = train_numeric,
    test_numeric  = test_numeric,
    train_factor  = train_factor,
    test_factor   = test_factor,
    train_lightgbm = train_lightgbm,
    test_lightgbm  = test_lightgbm,
    x_train_numeric = x_train_numeric,
    x_test_numeric  = x_test_numeric,
    train_indices   = train_indices,
    processed_data  = processed_data,
    lightgbm_model  = models_list$lightgbm
  )
}

#===================== 5. Counterfactual prediction function (prediction difference) =====================#

# -------- Fix 3: BAA_Reduction = predicted current - predicted counterfactual -------- #
counterfactual_mixed <- function(model_results, processed_data,
                                 target_evi = 0.40, target_o3 = 60) {
  original_data  <- processed_data$numeric_data
  best_model_name <- model_results$best_model_name
  cat("Performing counterfactual prediction with best model:", best_model_name, "\n")
  
  # 1. Construct counterfactual data (numeric)
  counterfactual_numeric <- original_data
  
  evi_need_intervention <- original_data$EVI_1KM_mean4y < target_evi
  o3_need_intervention  <- original_data$O3_1KM_mean4y  > target_o3
  
  counterfactual_numeric$EVI_1KM_mean4y[evi_need_intervention] <- target_evi
  counterfactual_numeric$O3_1KM_mean4y[o3_need_intervention]   <- target_o3
  
  if (best_model_name == "GBM") {
    counterfactual_factor <- processed_data$factor_data
    counterfactual_factor$EVI_1KM_mean4y[evi_need_intervention] <- target_evi
    counterfactual_factor$O3_1KM_mean4y[o3_need_intervention]   <- target_o3
  }
  if (best_model_name == "LightGBM") {
    counterfactual_lightgbm <- processed_data$lightgbm_data
    counterfactual_lightgbm$EVI_1KM_mean4y[evi_need_intervention] <- target_evi
    counterfactual_lightgbm$O3_1KM_mean4y[o3_need_intervention]   <- target_o3
  }
  
  best_model <- model_results$models[[tolower(best_model_name)]]
  
  # 2. Predict "current scenario" and "counterfactual scenario"
  if (best_model_name %in% c("LASSO", "ElasticNet")) {
    X_orig <- as.matrix(original_data %>% dplyr::select(-BAA_PhenoAge))
    X_cf   <- as.matrix(counterfactual_numeric %>% dplyr::select(-BAA_PhenoAge))
    
    pred_orig <- predict(best_model, newx = X_orig, s = "lambda.min")
    pred_cf   <- predict(best_model, newx = X_cf,   s = "lambda.min")
    
  } else if (best_model_name == "GBM") {
    best_trees <- gbm.perf(best_model, method = "cv", plot.it = FALSE)
    pred_orig  <- predict(best_model, processed_data$factor_data, n.trees = best_trees)
    pred_cf    <- predict(best_model, counterfactual_factor,       n.trees = best_trees)
    
  } else if (best_model_name == "LightGBM") {
    X_orig <- as.matrix(processed_data$lightgbm_data %>% dplyr::select(-BAA_PhenoAge))
    X_cf   <- as.matrix(counterfactual_lightgbm %>% dplyr::select(-BAA_PhenoAge))
    
    pred_orig <- predict(best_model, X_orig)
    pred_cf   <- predict(best_model, X_cf)
    
  } else if (best_model_name == "GAM") {
    pred_orig <- predict(best_model, original_data)
    pred_cf   <- predict(best_model, counterfactual_numeric)
    
  } else if (best_model_name %in% c("Cubist", "PLS")) {
    X_orig <- as.matrix(original_data %>% dplyr::select(-BAA_PhenoAge))
    X_cf   <- as.matrix(counterfactual_numeric %>% dplyr::select(-BAA_PhenoAge))
    
    pred_orig <- predict(best_model, X_orig)
    pred_cf   <- predict(best_model, X_cf)
    
  } else {
    pred_orig <- predict(best_model, original_data)
    pred_cf   <- predict(best_model, counterfactual_numeric)
  }
  
  pred_orig <- as.numeric(pred_orig)
  pred_cf   <- as.numeric(pred_cf)
  
  # 3. Counterfactual improvement = predicted current - predicted counterfactual
  baa_reduction <- pred_orig - pred_cf
  
  results <- data.frame(
    ID = 1:nrow(original_data),
    Observed_BAA        = original_data$BAA_PhenoAge,
    Pred_BAA_current    = pred_orig,
    Pred_BAA_counterf   = pred_cf,
    BAA_Reduction       = baa_reduction,
    Original_EVI = original_data$EVI_1KM_mean4y,
    Original_O3  = original_data$O3_1KM_mean4y,
    age = original_data$age,
    evi_intervention = evi_need_intervention,
    o3_intervention  = o3_need_intervention
  )
  
  list(
    results = results,
    best_model_name = best_model_name
  )
}

#===================== 6. Main analysis function =====================#

main_mixed_analysis <- function(data) {
  cat("=== Mixed data type ML counterfactual prediction analysis ===\n\n")
  
  cat("1. Data preprocessing...\n")
  processed_data <- preprocess_ml_data_mixed(data)
  
  cat("2. Model training...\n")
  model_results <- train_ml_models_mixed(processed_data)
  
  cat("\n3. Counterfactual prediction...\n")
  cat("Counterfactual scenario: EVI raised to 0.40, O3 reduced to 60 μg/m³\n\n")
  counterfactual_results <- counterfactual_mixed(model_results, processed_data)
  
  # Summary
  results_summary <- counterfactual_results$results %>%
    summarise(
      N = n(),
      # Compare using means of predicted BAA
      Mean_Original_BAA      = mean(Pred_BAA_current),
      Mean_Counterfactual_BAA = mean(Pred_BAA_counterf),
      Mean_Reduction = mean(BAA_Reduction),
      SD_Reduction   = sd(BAA_Reduction),
      Median_Reduction = median(BAA_Reduction),
      Min_Reduction     = min(BAA_Reduction),
      Max_Reduction     = max(BAA_Reduction),
      P25_Reduction     = quantile(BAA_Reduction, 0.25),
      P75_Reduction     = quantile(BAA_Reduction, 0.75),
      Proportion_Positive = mean(BAA_Reduction > 0),
      CI_Lower = Mean_Reduction - 1.96 * (SD_Reduction / sqrt(N)),
      CI_Upper = Mean_Reduction + 1.96 * (SD_Reduction / sqrt(N))
    )
  cat("4. Counterfactual prediction results summary:\n")
  print(results_summary)
  
  age_stratified <- counterfactual_results$results %>%
    mutate(Age_Group = cut(
      age,
      breaks = c(0, 60, 70, 80, 100),
      labels = c("<=60 years", "61-70 years", "71-80 years", ">80 years")
    )) %>%
    group_by(Age_Group) %>%
    summarise(
      N = n(),
      Mean_Reduction = mean(BAA_Reduction),
      SE_Reduction   = sd(BAA_Reduction) / sqrt(n()),
      CI_Lower = Mean_Reduction - 1.96 * SE_Reduction,
      CI_Upper = Mean_Reduction + 1.96 * SE_Reduction,
      .groups = "drop"
    )
  cat("\n5. Age-stratified analysis:\n")
  print(age_stratified)
  
  intervention_stats <- counterfactual_results$results %>%
    summarise(
      evi_intervention_pct   = mean(evi_intervention) * 100,
      o3_intervention_pct    = mean(o3_intervention) * 100,
      both_intervention_pct   = mean(evi_intervention & o3_intervention) * 100,
      any_intervention_pct   = mean(evi_intervention | o3_intervention) * 100
    )
  cat("\n6. Intervention statistics:\n")
  print(intervention_stats)
  
  list(
    processed_data = processed_data,
    model_results = model_results,
    counterfactual_results = counterfactual_results,
    results_summary = results_summary,
    age_stratified = age_stratified,
    intervention_stats = intervention_stats
  )
}

#===================== 7. Visualization functions (original versions retained) =====================#

# ---- 7.1 Comprehensive model performance visualization ----
create_comprehensive_performance_plots <- function(performance_df) {
  plots <- list()
  
  # 1. Main metric comparison (R2, RMSE, MAE)
  main_metrics <- performance_df %>%
    dplyr::select(Model, R2, RMSE, MAE) %>%
    tidyr::pivot_longer(cols = c(R2, RMSE, MAE),
                        names_to = "Metric", values_to = "Value")
  
  scale_metrics <- main_metrics %>%
    group_by(Metric) %>%
    mutate(Scaled_Value = ifelse(
      Metric == "R2",
      Value,
      (max(Value) - Value) / (max(Value) - min(Value))
    ))
  
  p1 <- ggplot(scale_metrics, aes(x = Model, y = Scaled_Value, fill = Metric)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_text(aes(label = round(Value, 3)),
              position = position_dodge(0.9), vjust = -0.5, size = 2.5) +
    scale_fill_manual(values = c("R2" = "#2E8B57", "RMSE" = "#FF6347", "MAE" = "#1E90FF")) +
    labs(title = "Comparison of main performance metrics",
         subtitle = "R2 (higher is better); RMSE and MAE (standardized, higher is better)",
         y = "Standardized metric value", x = "Model") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plots$main_metrics <- p1
  
  # 2. R² vs Adj_R2
  r2_data <- performance_df %>%
    dplyr::select(Model, R2, Adj_R2) %>%
    tidyr::pivot_longer(cols = c(R2, Adj_R2),
                        names_to = "R2_Type", values_to = "Value")
  p2 <- ggplot(r2_data, aes(x = Model, y = Value, fill = R2_Type)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_text(aes(label = round(Value, 4)),
              position = position_dodge(0.9), vjust = -0.5, size = 2.5) +
    scale_fill_manual(values = c("R2" = "#2E8B57", "Adj_R2" = "#32CD32")) +
    labs(title = "R2 and adjusted R2 comparison", y = "R2 value", x = "Model") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plots$r2_comparison <- p2
  
  # 3. Performance ranking heatmap
  rank_data <- performance_df %>%
    dplyr::select(Model, R2, RMSE, MAE, Explained_Variance, MAPE) %>%
    mutate(
      R2_Rank  = rank(-R2),
      RMSE_Rank = rank(RMSE),
      MAE_Rank  = rank(MAE),
      EV_Rank   = rank(-Explained_Variance),
      MAPE_Rank = rank(MAPE)
    ) %>%
    dplyr::select(Model, R2_Rank, RMSE_Rank, MAE_Rank, EV_Rank, MAPE_Rank) %>%
    tidyr::pivot_longer(cols = -Model, names_to = "Metric", values_to = "Rank")
  p3 <- ggplot(rank_data, aes(x = Metric, y = Model, fill = Rank)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = Rank), color = "white", size = 4, fontface = "bold") +
    scale_fill_gradient(low = "#DC143C", high = "#32CD32", name = "Rank\n(1=best)") +
    labs(title = "Model performance ranking heatmap",
         subtitle = "Ranking across multiple metrics (1 = best)",
         x = "Evaluation metric", y = "Model") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plots$rank_heatmap <- p3
  
  # 4. Composite scores
  score_data <- performance_df %>%
    mutate(
      R2_Score  = (R2  - min(R2))  / (max(R2)  - min(R2)),
      RMSE_Score = 1 - (RMSE - min(RMSE))/(max(RMSE) - min(RMSE)),
      MAE_Score  = 1 - (MAE  - min(MAE)) /(max(MAE)  - min(MAE)),
      EV_Score   = (Explained_Variance - min(Explained_Variance)) /
        (max(Explained_Variance) - min(Explained_Variance)),
      Composite_Score = (R2_Score + RMSE_Score + MAE_Score + EV_Score) / 4
    ) %>%
    dplyr::select(Model, Composite_Score, R2_Score, RMSE_Score, MAE_Score, EV_Score) %>%
    tidyr::pivot_longer(cols = -Model, names_to = "Score_Type", values_to = "Score")
  p4 <- ggplot(score_data, aes(x = Model, y = Score, fill = Score_Type)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_text(aes(label = round(Score, 3)),
              position = position_dodge(0.9), vjust = -0.5, size = 2.5) +
    scale_fill_manual(values = c(
      "Composite_Score" = "#FF4500",
      "R2_Score"  = "#2E8B57",
      "RMSE_Score"= "#1E90FF",
      "MAE_Score" = "#9370DB",
      "EV_Score"  = "#FFD700"
    )) +
    labs(title = "Composite performance scores",
         subtitle = "Normalized scores across metrics (0-1, 1 = best)",
         y = "Normalized score", x = "Model") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plots$composite_scores <- p4
  
  # 5. Parallel coordinates plot
  parallel_data <- performance_df %>%
    dplyr::select(Model, R2, RMSE, MAE, Explained_Variance) %>%
    mutate(across(c(R2, Explained_Variance),
                  ~ (. - min(.)) / (max(.) - min(.)),
                  .names = "Scaled_{.col}")) %>%
    mutate(across(c(RMSE, MAE),
                  ~ 1 - (. - min(.)) / (max(.) - min(.)),
                  .names = "Scaled_{.col}"))
  
  parallel_long <- parallel_data %>%
    dplyr::select(Model, Scaled_R2, Scaled_RMSE, Scaled_MAE, Scaled_Explained_Variance) %>%
    tidyr::pivot_longer(cols = -Model, names_to = "Metric", values_to = "Scaled_Value")
  
  p5 <- ggplot(parallel_long,
               aes(x = Metric, y = Scaled_Value,
                   color = Model, group = Model)) +
    geom_line(size = 1.5, alpha = 0.7) +
    geom_point(size = 3) +
    labs(title = "Parallel coordinates plot of model performance",
         subtitle = "All metrics normalized to 0-1 (1 = best)",
         y = "Standardized metric value", x = "Evaluation metric") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  plots$parallel_coordinates <- p5
  
  plots
}

# ---- 7.2 Main results figures (BAA improvement distribution etc.) ----
# ---- Overwrite the original create_visualizations function ----
create_visualizations <- function(analysis_results) {
  
  # Common large-font academic theme
  base_theme <- theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(
        family = "serif", face = "bold",
        size = 20, hjust = 0.5, margin = margin(b = 6)
      ),
      plot.subtitle = element_text(
        size = 12, hjust = 0.5,
        margin = margin(b = 3)
      ),
      plot.caption = element_text(
        size = 10, hjust = 0.5,
        margin = margin(t = 6)
      ),
      axis.title.x = element_text(
        size = 16, face = "bold",
        margin = margin(t = 6)
      ),
      axis.title.y = element_text(
        size = 16, face = "bold",
        margin = margin(r = 8)
      ),
      axis.text.x  = element_text(size = 13),
      axis.text.y  = element_text(size = 13),
      panel.grid.major.y = element_line(
        colour = "grey90", linetype = "dashed"
      ),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 12),
      plot.margin = margin(6, 10, 6, 6)
    )
  
  dat <- analysis_results$counterfactual_results$results
  
  ## 1. Distribution of BAA improvement ----
  # Calculate max histogram height to reserve space for text annotations
  hobj   <- hist(dat$BAA_Reduction, breaks = 40, plot = FALSE)
  y_max  <- max(hobj$counts, na.rm = TRUE)
  y_text <- y_max * 1.18
  
  p1 <- ggplot(dat, aes(x = BAA_Reduction)) +
    geom_histogram(
      bins = 40,
      fill = "steelblue", alpha = 0.75,
      color = "white"
    ) +
    geom_vline(
      xintercept = 0,
      color = "red", linetype = "dashed", linewidth = 1
    ) +
    geom_vline(
      xintercept = analysis_results$results_summary$Mean_Reduction,
      color = "navy", linetype = "dashed", linewidth = 1
    ) +
    annotate(
      "text",
      x = analysis_results$results_summary$Mean_Reduction,
      y = y_text,
      label = paste0("Mean improvement: ",
                     round(analysis_results$results_summary$Mean_Reduction, 2),
                     " years"),
      vjust = 0, size = 4.5, colour = "navy", fontface = "bold"
    ) +
    coord_cartesian(ylim = c(0, y_max * 1.30)) +
    labs(
      title    = "Distribution of BAA improvement under the counterfactual scenario",
      subtitle = paste0("Best-performing model: ",
                        analysis_results$counterfactual_results$best_model_name),
      x        = "BAA improvement (years)",
      y        = "Number of participants",
      caption  = "Positive values indicate improvement (lower biological age)."
    ) +
    base_theme
  
  ## 2. Mean BAA improvement by age group ----
  age_data <- analysis_results$age_stratified
  age_data <- age_data[complete.cases(age_data), ]
  
  p2 <- ggplot(age_data, aes(x = Age_Group, y = Mean_Reduction)) +
    geom_col(fill = "#F8766D", alpha = 0.8, width = 0.7) +
    geom_errorbar(
      aes(ymin = CI_Lower, ymax = CI_Upper),
      width = 0.15, linewidth = 0.7
    ) +
    # Place labels above error bar tops
    geom_text(
      aes(y = CI_Upper, label = round(Mean_Reduction, 2)),
      vjust = -0.8, size = 4.2
    ) +
    labs(
      title = "Mean BAA improvement by age group",
      x     = "Age group",
      y     = "Mean BAA improvement (years)"
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(size = 14, angle = 0, vjust = 0.5)
    ) +
    coord_cartesian(ylim = c(
      min(age_data$CI_Lower, na.rm = TRUE) * 0.9,
      max(age_data$CI_Upper, na.rm = TRUE) * 1.20
    ))
  
  ## 3. Performance plots (tweak style only, no data changes) ----
  performance_plots <- create_comprehensive_performance_plots(
    analysis_results$model_results$performance
  )
  
  performance_plots$main_metrics <- performance_plots$main_metrics +
    labs(
      title    = "Comparison of main performance metrics",
      subtitle = "R² (higher is better); RMSE and MAE (standardized, higher is better)",
      x        = "Model",
      y        = "Standardized metric value"
    ) + base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  performance_plots$r2_comparison <- performance_plots$r2_comparison +
    labs(
      title = "R² and adjusted R² across models",
      x     = "Model",
      y     = "Coefficient of determination"
    ) + base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  performance_plots$rank_heatmap <- performance_plots$rank_heatmap +
    labs(
      title    = "Ranking of model performance across metrics",
      subtitle = "Rank 1 indicates the best performance",
      x        = "Evaluation metric",
      y        = "Model"
    ) + base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  performance_plots$composite_scores <- performance_plots$composite_scores +
    labs(
      title    = "Composite performance scores of models",
      subtitle = "All metrics normalized to a 0–1 scale (1 = best)",
      x        = "Model",
      y        = "Normalized score"
    ) + base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  performance_plots$parallel_coordinates <- performance_plots$parallel_coordinates +
    labs(
      title    = "Parallel coordinate plot of model performance",
      subtitle = "All metrics normalized to a 0–1 scale (1 = best)",
      x        = "Evaluation metric",
      y        = "Normalized metric value"
    ) + base_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  ## 4. Baseline EVI vs BAA improvement ----
  max_imp <- max(dat$BAA_Reduction, na.rm = TRUE)
  
  p4 <- ggplot(dat, aes(x = Original_EVI, y = BAA_Reduction)) +
    geom_point(alpha = 0.6, color = "darkgreen", size = 1.6) +
    geom_smooth(method = "lm", color = "red", linewidth = 0.8, se = TRUE) +
    geom_vline(xintercept = 0.40, linetype = "dashed", color = "navy") +
    annotate(
      "text",
      x = 0.40, y = max_imp * 0.95,
      label = "Target EVI = 0.40",
      hjust = -0.05, size = 4.2, colour = "navy", fontface = "bold"
    ) +
    labs(
      title = "Association of baseline greenness with BAA improvement",
      x     = "Baseline EVI",
      y     = "BAA improvement (years)"
    ) +
    base_theme
  
  ## 5. Baseline O3 vs BAA improvement ----
  p5 <- ggplot(dat, aes(x = Original_O3, y = BAA_Reduction)) +
    geom_point(alpha = 0.6, color = "#E69F00", size = 1.6) +
    geom_smooth(method = "lm", color = "red", linewidth = 0.8, se = TRUE) +
    geom_vline(xintercept = 60, linetype = "dashed", color = "navy") +
    annotate(
      "text",
      x = 60, y = max_imp * 0.95,
      label = "Target O3 = 60 µg/m³",
      hjust = -0.05, size = 4.2, colour = "navy", fontface = "bold"
    ) +
    labs(
      title = "Association of baseline ozone with BAA improvement",
      x     = "Baseline O3 concentration (µg/m³)",
      y     = "BAA improvement (years)"
    ) +
    base_theme
  
  ## 6. Single-population boxplot (optional, currently unused in combined figure) ----
  p6 <- ggplot(dat, aes(x = "", y = BAA_Reduction)) +
    geom_boxplot(
      fill = "#00BFC4", alpha = 0.75,
      width = 0.35, outlier.shape = NA
    ) +
    stat_summary(
      fun = mean, geom = "point",
      shape = 18, size = 3.5, color = "red"
    ) +
    labs(
      title = "Distribution of individual-level BAA improvement",
      x     = "",
      y     = "BAA improvement (years)"
    ) +
    base_theme +
    theme(axis.text.x = element_blank())
  
  ## 7. Density plot (backup) ----
  p7 <- ggplot(dat, aes(x = BAA_Reduction)) +
    geom_density(fill = "lightgreen", alpha = 0.7) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    geom_vline(
      xintercept = analysis_results$results_summary$Mean_Reduction,
      color = "navy", linetype = "dashed"
    ) +
    labs(
      title = "Density distribution of BAA improvement",
      x     = "BAA improvement (years)",
      y     = "Density"
    ) +
    base_theme
  
  list(
    improvement_dist     = p1,
    age_stratified       = p2,
    performance_plots    = performance_plots,
    evi_effect           = p4,
    o3_effect            = p5,
    improvement_boxplot  = p6,
    density_plot         = p7
  )
}


#===================== 8. LightGBM SHAP functions =====================#
# These functions follow the previous version exactly (create_lightgbm_shap_analysis_full,
# create_shap_explainability_plot, create_combined_importance_shap_plot,
# create_advanced_shap_summary_plot) - kept verbatim without additional explanation,
# only assembled together.

# ---- 8.1 LightGBM full-sample SHAP ----
create_lightgbm_shap_analysis_full <- function(model_results) {
  cat("Performing LightGBM model SHAP analysis (full sample)...\n")
  
  lightgbm_model <- model_results$lightgbm_model
  if (is.null(lightgbm_model)) {
    cat("LightGBM model not found, skipping SHAP analysis\n")
    return(NULL)
  }
  
  test_data <- model_results$test_lightgbm
  x_test <- test_data %>% dplyr::select(-BAA_PhenoAge)
  y_test <- test_data$BAA_PhenoAge
  
  cat("Test data dimensions:", dim(x_test), "\n")
  cat("Variable names:", names(x_test), "\n")
  cat("Full-sample SHAP analysis, sample size:", nrow(x_test), "\n")
  
  x_test_sample <- as.matrix(x_test)
  y_test_sample <- y_test
  
  tryCatch({
    pred_wrapper <- function(object, newdata) {
      predict(object, newdata)
    }
    
    baseline <- mean(predict(lightgbm_model, x_test_sample))
    cat("Baseline prediction:", baseline, "\n")
    
    cat("Starting SHAP value computation (full sample, nsim = 10)...\n")
    shap_values <- fastshap::explain(
      lightgbm_model,
      X = x_test_sample,
      pred_wrapper = pred_wrapper,
      nsim = 10
    )
    cat("SHAP value computation complete, dimensions:", dim(shap_values), "\n")
    
    # 1. SHAP summary plot
    cat("Generating SHAP summary plot...\n")
    shap_importance <- colMeans(abs(shap_values)) %>%
      sort(decreasing = TRUE) %>%
      head(15)
    
    plot_data <- data.frame()
    for (feature in names(shap_importance)) {
      temp_df <- data.frame(
        feature = feature,
        shap_value = shap_values[, feature],
        feature_value = x_test_sample[, feature]
      )
      plot_data <- rbind(plot_data, temp_df)
    }
    
    p_shap_summary <- ggplot(plot_data,
                             aes(x = shap_value,
                                 y = reorder(feature, abs(shap_value), mean))) +
      geom_point(aes(color = feature_value), alpha = 0.6, size = 1.2) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_viridis_c(name = "Feature value") +
      labs(
        title = "LightGBM model SHAP summary plot",
        subtitle = paste("Based on", nrow(x_test), "samples' SHAP analysis"),
        x = "SHAP value (impact on prediction)",
        y = "Feature"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "right"
      )
    
    # 2. SHAP waterfall plot (randomly select 4 samples)
    cat("Generating SHAP waterfall plot...\n")
    set.seed(20240102)
    waterfall_samples <- sample(1:nrow(x_test_sample),
                                min(4, nrow(x_test_sample)))
    waterfall_plots <- list()
    
    for (i in seq_along(waterfall_samples)) {
      idx <- waterfall_samples[i]
      sample_shap <- shap_values[idx, ]
      sample_features <- x_test_sample[idx, ]
      pred_value <- predict(lightgbm_model,
                            x_test_sample[idx, , drop = FALSE])
      
      waterfall_data <- data.frame(
        feature = names(sample_shap),
        shap_value = as.numeric(sample_shap)
      ) %>%
        mutate(abs_shap = abs(shap_value)) %>%
        arrange(desc(abs_shap)) %>%
        head(8)
      
      p_waterfall <- ggplot(waterfall_data,
                            aes(x = reorder(feature, abs_shap),
                                y = shap_value)) +
        geom_col(aes(fill = shap_value > 0), alpha = 0.8) +
        coord_flip() +
        scale_fill_manual(values = c("TRUE" = "#FF6B6B",
                                     "FALSE" = "#4ECDC4"),
                          labels = c("TRUE" = "Positive", "FALSE" = "Negative"),
                          name   = "Direction") +
        labs(
          title = paste("SHAP waterfall - Sample", idx),
          subtitle = paste("Predicted: ", round(pred_value, 2),
                           "Actual: ", round(y_test_sample[idx], 2)),
          x = "Feature", y = "SHAP value contribution"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 12),
          axis.text.y = element_text(size = 9)
        )
      waterfall_plots[[i]] <- p_waterfall
    }
    
    p_waterfall_combined <- wrap_plots(waterfall_plots, ncol = 2) +
      plot_annotation(
        title = "LightGBM model SHAP waterfall examples",
        theme = theme(plot.title = element_text(face = "bold", size = 14))
      )
    
    list(
      shap_summary = p_shap_summary,
      waterfall    = p_waterfall_combined,
      shap_values  = shap_values,
      baseline     = baseline
    )
  }, error = function(e) {
    cat("SHAP analysis failed:", e$message, "\n")
    NULL
  })
}

# ---- 8.2 SHAP + density-weighted beeswarm ----
# ---- 8.2 SHAP + density-weighted beeswarm (EN, large fonts + EN y-axis labels) ----
create_shap_explainability_plot <- function(shap_results, model_results, output_dir) {
  cat("Creating SHAP explainability plot (beeswarm + density weighting)...\n")
  if (is.null(shap_results)) {
    cat("SHAP results are NULL, skip.\n"); return(NULL)
  }
  
  test_data <- model_results$test_lightgbm
  if (is.null(test_data)) {
    cat("Test data not found, skip.\n"); return(NULL)
  }
  x_test_for_plot <- test_data %>% dplyr::select(-BAA_PhenoAge)
  
  # ---------- helper: robust rescaling 0–1 ----------
  robust01 <- function(x, lo = 0.05, hi = 0.95) {
    ql <- quantile(x, probs = lo, na.rm = TRUE)
    qh <- quantile(x, probs = hi, na.rm = TRUE)
    if (qh <= ql) return(rep(0.5, length(x)))
    pmin(1, pmax(0, (x - ql) / (qh - ql)))
  }
  
  # ---------- core drawing function ----------
  make_shap_combined_plot_ref <- function(shap_values, X, top_k = 15,
                                          title_top = "Mean(|SHAP|) across participants",
                                          align_at = -4) {
    S <- as.matrix(shap_values)
    stopifnot(nrow(S) == nrow(X))
    
    feat_names <- colnames(S)
    if (is.null(feat_names)) feat_names <- paste0("f", seq_len(ncol(S)))
    colnames(S) <- feat_names
    
    # rank by mean(|SHAP|)
    mabs <- colMeans(abs(S), na.rm = TRUE)
    ord  <- order(mabs, decreasing = TRUE)
    feat_sorted <- feat_names[ord]
    mabs_sorted <- mabs[ord]
    if (length(feat_sorted) > top_k) {
      feat_sorted <- feat_sorted[1:top_k]
      mabs_sorted <- mabs_sorted[1:top_k]
    }
    pct_sorted <- mabs_sorted / sum(mabs_sorted) * 100
    
    # ---------- map original variable names to English labels ----------
    var_label_map <- c(
      "EVI_1KM_mean4y" = "EVI",
      "O3_1KM_mean4y"  = "O\u2083",                  # O₃ (Unicode subscript 3)
      "age"            = "Age (years)",
      "sex"            = "Sex",
      "ln_hh_consumption_total"     = "Log household consumption",
      "retirement"         = "Retired",
      "pension_insurance"         = "Pension insurance",
      "smoking_current"         = "Current smoking",
      "alcohol_current"         = "Current alcohol drinking",
      "climate_pc1"     = "Climate PC1",
      "climate_pc2"     = "Climate PC2",
      "climate_pc3"     = "Climate PC3",
      "climate_pc4"     = "Climate PC4",
      "climate_pc5"     = "Climate PC5",
      "bmi"      = "Body mass index",
      "hypertension"           = "Hypertension",
      "dyslipidemia"         = "Dyslipidemia"
    )
    
    english_labels_sorted <- vapply(
      feat_sorted,
      function(f) if (!is.na(var_label_map[f])) var_label_map[f] else f,
      character(1)
    )
    
    Xsel <- as.data.frame(X)[, feat_sorted, drop = FALSE]
    Ssel <- S[, feat_sorted, drop = FALSE]
    
    # beeswarm points with density-weighted jitter
    bees_list <- vector("list", length(feat_sorted))
    for (i in seq_along(feat_sorted)) {
      fv  <- Xsel[[i]]
      col <- robust01(fv)
      sv  <- Ssel[, i]
      
      dens <- density(sv[is.finite(sv)], n = 512, na.rm = TRUE)
      dval <- approx(dens$x, dens$y, xout = sv, rule = 2, ties = "ordered")$y
      dval[!is.finite(dval)] <- 0
      dscaled <- if (max(dval) > 0) dval / max(dval) else dval
      yoff <- rnorm(length(sv), sd = 0.15) * dscaled
      
      bees_list[[i]] <- data.frame(
        feature = feat_sorted[i],
        shap_value = sv,
        feature_value_col = col,
        y_offset = yoff
      )
    }
    
    bees_long <- dplyr::bind_rows(bees_list)
    bees_long$feature <- factor(bees_long$feature, levels = rev(feat_sorted))
    base_y <- as.numeric(bees_long$feature)
    bees_long$y <- base_y + pmin(0.48, pmax(-0.48, bees_long$y_offset))
    
    # horizontal bars for mean(|SHAP|)
    bars_df <- data.frame(
      feature = factor(feat_sorted, levels = rev(feat_sorted)),
      mabs    = mabs_sorted,
      pct     = pct_sorted,
      xmin    = align_at
    )
    max_mabs <- max(bars_df$mabs, na.rm = TRUE)
    k <- if (max_mabs == 0) 1 else max(abs(bees_long$shap_value), na.rm = TRUE) / max_mabs
    bars_df$xmax <- bars_df$xmin + bars_df$mabs * k
    
    label_df <- bars_df %>%
      dplyr::mutate(
        label = sprintf("%.3f (%.1f%%)", mabs, pct),
        x = xmin + 0.02 * k * max_mabs,
        y = as.numeric(feature)
      )
    
    col_low  <- "#3b4cc0"
    col_high <- "#b40486"
    
    right_edge_points <- max(bees_long$shap_value, na.rm = TRUE)
    right_edge_bars   <- max(bars_df$xmax, na.rm = TRUE)
    x_right <- max(right_edge_points, right_edge_bars)
    xrange  <- x_right - align_at
    xlim_c  <- c(align_at, align_at + 1.05 * xrange)
    
    # scientific theme (large fonts)
    theme_scientific <- theme_classic(base_size = 18) +
      theme(
        plot.title      = element_text(face = "bold", size = 24, hjust = 0.5),
        axis.title.x    = element_text(size = 20, face = "bold",
                                       margin = margin(t = 10)),
        axis.title.y    = element_text(size = 20, face = "bold",
                                       margin = margin(r = 10)),
        axis.text.x     = element_text(size = 16),
        axis.text.y     = element_text(size = 16),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "right",
        legend.title       = element_text(size = 18, face = "bold"),
        legend.text        = element_text(size = 16),
        axis.title.x.top   = element_text(size = 18, face = "bold",
                                          margin = margin(b = 8))
      )
    
    p <- ggplot() +
      # bars for mean(|SHAP|)
      geom_rect(
        data = bars_df,
        aes(ymin = as.numeric(feature) - 0.48,
            ymax = as.numeric(feature) + 0.48,
            xmin = xmin, xmax = xmax),
        fill = "#f6b6c5", alpha = 0.45, color = NA
      ) +
      # labels for mean(|SHAP|)
      geom_label(
        data = label_df,
        aes(x = x, y = y, label = label),
        linewidth = 0,
        label.padding = unit(0.10, "lines"),
        label.r = unit(0.06, "lines"),
        fill = alpha("white", 0.9),
        hjust = 0, vjust = 0.5, size = 4.2, colour = "#222222"
      ) +
      # beeswarm points
      geom_point(
        data = bees_long,
        aes(x = shap_value, y = y, colour = feature_value_col),
        size = 1.8, alpha = 0.75
      ) +
      geom_vline(xintercept = 0, colour = "grey40", linetype = "dashed") +
      scale_color_gradient(
        low   = col_low,
        high  = col_high,
        name  = "Feature value",
        limits = c(0, 1),
        breaks = c(0, 1),
        labels = c("Low", "High"),
        guide = guide_colorbar(
          title.position = "top",
          barheight = unit(18, "lines")
        )
      ) +
      scale_x_continuous(
        name   = "SHAP value (impact on model prediction)",
        limits = xlim_c,
        sec.axis = sec_axis(
          ~ (. - align_at) / k,
          breaks = pretty(c(0, max_mabs)),
          labels = scales::number_format(accuracy = 0.02),
          name   = title_top
        )
      ) +
      scale_y_continuous(
        name   = NULL,
        breaks = seq_along(rev(feat_sorted)),
        labels = rev(english_labels_sorted)   # <── English feature names + O₃
      ) +
      ggtitle("Mean absolute SHAP values and individual contributions") +
      theme_scientific
    
    p
  }
  
  # ---------- build & save figure ----------
  tryCatch({
    p_combined <- make_shap_combined_plot_ref(
      shap_values = shap_results$shap_values,
      X           = x_test_for_plot,
      top_k       = 15,
      title_top   = "Mean(|SHAP|) across participants",
      align_at    = -4
    ) +
      labs(
        title    = "SHAP explainability plot (density-weighted beeswarm)",
        subtitle = "Based on the LightGBM model for BAA prediction"
      )
    
    shap_dir <- file.path(output_dir, "SHAP_analysis")
    if (!dir.exists(shap_dir)) dir.create(shap_dir, recursive = TRUE)
    outfile <- file.path(shap_dir,
                         "SHAP_explainability_beeswarm.png")
    
    ggsave(outfile, p_combined, width = 14, height = 7.5, dpi = 300)
    cat("✓ SHAP explainability plot saved to:", outfile, "\n")
    print(p_combined)
    p_combined
  }, error = function(e) {
    cat("✗ Failed to create SHAP explainability plot:", e$message, "\n")
    NULL
  })
}



# ---- 8.3 LightGBM feature importance vs SHAP ----
create_combined_importance_shap_plot <- function(shap_results, model_results, output_dir) {
  cat("Creating feature importance vs SHAP combined plot...\n")
  if (is.null(shap_results)) {
    cat("SHAP results are empty\n"); return(NULL)
  }
  
  lightgbm_model <- model_results$lightgbm_model
  if (is.null(lightgbm_model)) {
    cat("LightGBM model not found\n"); return(NULL)
  }
  
  tryCatch({
    importance_matrix <- lgb.importance(lightgbm_model)
    if (is.null(importance_matrix) || nrow(importance_matrix) == 0) {
      cat("Cannot retrieve feature importance\n"); return(NULL)
    }
    
    shap_values <- shap_results$shap_values
    mean_abs_shap <- colMeans(abs(shap_values))
    shap_importance <- data.frame(
      Feature = names(mean_abs_shap),
      mean_SHAP = mean_abs_shap,
      stringsAsFactors = FALSE
    )
    
    combined_importance <- importance_matrix %>%
      dplyr::select(Feature, Gain) %>%
      left_join(shap_importance, by = "Feature") %>%
      mutate(
        Gain_norm = Gain / max(Gain, na.rm = TRUE),
        SHAP_norm = mean_SHAP / max(mean_SHAP, na.rm = TRUE)
      ) %>%
      arrange(desc(Gain_norm + SHAP_norm)) %>%
      head(15)
    
    combined_long <- combined_importance %>%
      dplyr::select(Feature, Gain_norm, SHAP_norm) %>%
      tidyr::pivot_longer(
        cols = c(Gain_norm, SHAP_norm),
        names_to = "Importance_Type",
        values_to = "Normalized_Importance"
      ) %>%
      mutate(
        Importance_Type = factor(
          Importance_Type,
          levels = c("Gain_norm", "SHAP_norm"),
          labels = c("Feature importance (Gain)", "SHAP importance")
        )
      )
    
    p_combined <- ggplot(combined_long,
                         aes(x = reorder(Feature, Normalized_Importance),
                             y = Normalized_Importance,
                             fill = Importance_Type)) +
      geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
      coord_flip() +
      scale_fill_manual(
        values = c("Feature importance (Gain)" = "#FF6B6B",
                   "SHAP importance"       = "#4ECDC4"),
        name = "Importance type"
      ) +
      labs(
        title = "LightGBM model: feature importance vs SHAP importance",
        subtitle = "Top 15 features (normalized to 0-1)",
        x = "Feature", y = "Normalized importance"
      ) +
      theme_minimal() +
      theme(
        plot.title  = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50"),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")
      ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1))
    
    shap_dir <- file.path(output_dir, "SHAP_analysis")
    if (!dir.exists(shap_dir)) dir.create(shap_dir, recursive = TRUE)
    ggsave(file.path(shap_dir, "LightGBM_importance_vs_SHAP.png"),
           p_combined, width = 14, height = 10, dpi = 300)
    cat("✓ Feature importance vs SHAP combined plot saved successfully\n")
    print(p_combined)
    p_combined
  }, error = function(e) {
    cat("✗ Feature importance vs SHAP combined plot failed:", e$message, "\n")
    NULL
  })
}

# ---- 8.4 Advanced SHAP summary plot ----
create_advanced_shap_summary_plot <- function(shap_results, model_results, output_dir) {
  cat("Creating advanced SHAP summary plot...\n")
  if (is.null(shap_results)) { cat("SHAP results are empty\n"); return(NULL) }
  
  tryCatch({
    shap_values <- shap_results$shap_values
    x_test_sample <- as.matrix(model_results$test_lightgbm %>% dplyr::select(-BAA_PhenoAge))
    if (is.null(shap_values) || is.null(x_test_sample)) {
      cat("SHAP values or test data not found\n"); return(NULL)
    }
    
    mean_abs_shap <- colMeans(abs(shap_values))
    feature_names <- names(mean_abs_shap)
    shap_df <- data.frame(feature = feature_names,
                          mean_shap = mean_abs_shap) %>%
      arrange(desc(mean_shap)) %>%
      head(15)
    sorted_features    <- shap_df$feature
    mean_shap_sorted   <- shap_df$mean_shap
    mean_shap_percent  <- mean_shap_sorted / sum(mean_shap_sorted) * 100
    
    plot_data <- data.frame()
    for (i in seq_along(sorted_features)) {
      feature <- sorted_features[i]
      temp_data <- data.frame(
        feature = feature,
        shap_value = shap_values[, feature],
        feature_value = x_test_sample[, feature],
        stringsAsFactors = FALSE
      )
      plot_data <- rbind(plot_data, temp_data)
    }
    plot_data$feature <- factor(plot_data$feature, levels = rev(sorted_features))
    
    p_advanced <- ggplot(plot_data, aes(x = shap_value, y = feature)) +
      geom_point(aes(color = feature_value), alpha = 0.6, size = 2, shape = 16) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
      geom_text(
        data = data.frame(
          feature = sorted_features,
          mean_shap = mean_shap_sorted,
          mean_percent = mean_shap_percent
        ),
        aes(x = -max(abs(plot_data$shap_value)) * 0.9,
            y = feature,
            label = sprintf("%.3f (%.1f%%)", mean_shap, mean_percent)),
        hjust = 0, size = 3.5, color = "black", fontface = "bold"
      ) +
      scale_color_viridis_c(
        name = "Feature value",
        option = "plasma",
        guide = guide_colorbar(barwidth = 12, barheight = 0.8,
                               direction = "horizontal", title.position = "top")
      ) +
      labs(
        title = "LightGBM model SHAP advanced summary plot",
        subtitle = "Points: individual SHAP values coloured by feature value; left: mean |SHAP| with percentage",
        x = "SHAP value (impact on prediction)",
        y = "Feature"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50",
                                     margin = margin(b = 15)),
        axis.title = element_text(face = "bold", size = 12),
        axis.text  = element_text(size = 10),
        axis.text.y = element_text(face = "bold", color = "black"),
        legend.position = "bottom",
        legend.title    = element_text(face = "bold", size = 10),
        legend.text     = element_text(size = 9),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        plot.margin      = margin(20, 20, 20, 20)
      ) +
      xlim(-max(abs(plot_data$shap_value)) * 1.1,
           max(abs(plot_data$shap_value)) * 1.1)
    
    shap_dir <- file.path(output_dir, "SHAP_analysis")
    if (!dir.exists(shap_dir)) dir.create(shap_dir, recursive = TRUE)
    ggsave(file.path(shap_dir, "LightGBM_SHAP_advanced_summary.png"),
           p_advanced, width = 16, height = 12, dpi = 300)
    cat("✓ Advanced SHAP summary plot saved successfully\n")
    print(p_advanced)
    p_advanced
  }, error = function(e) {
    cat("✗ Advanced SHAP summary plot failed:", e$message, "\n")
    NULL
  })
}

#===================== 9. Run main analysis =====================#

cat("Starting mixed data analysis with climate PCA...\n")
mixed_results <- main_mixed_analysis(data)

cat("\n=== In-depth model performance analysis ===\n")
best_model_name <- mixed_results$counterfactual_results$best_model_name
cat("Best model:", best_model_name, "\n")

cat("Generating visualizations...\n")
plots <- create_visualizations(mixed_results)

print(plots$improvement_dist)
print(plots$age_stratified)
print(plots$performance_plots$main_metrics)
print(plots$performance_plots$r2_comparison)
print(plots$performance_plots$rank_heatmap)
print(plots$performance_plots$composite_scores)
print(plots$performance_plots$parallel_coordinates)

cat("\n=== Main analysis results ===\n")
cat("Mean BAA improvement:", round(mixed_results$results_summary$Mean_Reduction, 3), " years\n")
cat("Proportion with improvement:", round(mixed_results$results_summary$Proportion_Positive * 100, 2), "%\n")
cat("Best model:", mixed_results$counterfactual_results$best_model_name,
    "(R²=",
    round(mixed_results$model_results$performance$R2[
      mixed_results$model_results$performance$Model ==
        mixed_results$counterfactual_results$best_model_name
    ], 4), ")\n")
cat("Proportion requiring environmental intervention:",
    round(mean(mixed_results$counterfactual_results$results$evi_intervention |
                 mixed_results$counterfactual_results$results$o3_intervention) * 100, 2), "%\n")

#===================== 10. Output directories (absolute paths) =====================#

project_root <- "."
output_fig_dir  <- "output/figures"
output_tab_dir  <- "output/tables"
sapply(c(output_fig_dir, output_tab_dir), function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE))

cat("Figure output:", output_fig_dir, "\n")
cat("Table output: ", output_tab_dir, "\n")

# ---- 10a. 4-panel counterfactual summary (intermediate component) ----
combined_plot_main <- (plots$improvement_dist + plots$age_stratified) /
  (plots$evi_effect + plots$o3_effect) +
  plot_layout(heights = c(1, 1))

# ---- 10b. Performance assessment plots ----
performance_combined <- (plots$performance_plots$main_metrics +
                           plots$performance_plots$r2_comparison) /
  (plots$performance_plots$rank_heatmap +
     plots$performance_plots$composite_scores) /
  plots$performance_plots$parallel_coordinates

tryCatch({
  ggsave(file.path(output_fig_dir, "model_performance_assessment_LightGBM.png"),
         performance_combined, width = 18, height = 16, dpi = 300)
  cat("Performance assessment combined plot saved successfully\n")
}, error = function(e) cat("Performance assessment combined plot save failed:", e$message, "\n"))

performance_charts <- list(
  "main_metrics_comparison_LightGBM.png"  = plots$performance_plots$main_metrics,
  "R2_vs_AdjR2_LightGBM.png"   = plots$performance_plots$r2_comparison,
  "model_performance_rank_heatmap_LightGBM.png" = plots$performance_plots$rank_heatmap,
  "model_composite_scores_LightGBM.png" = plots$performance_plots$composite_scores,
  "model_parallel_coordinates_LightGBM.png" = plots$performance_plots$parallel_coordinates
)

for (filename in names(performance_charts)) {
  tryCatch({
    ggsave(file.path(output_fig_dir, filename),
           performance_charts[[filename]], width = 12, height = 8, dpi = 300)
    cat("✓", filename, "saved successfully\n")
  }, error = function(e) cat("✗", filename, "save failed:", e$message, "\n"))
}

# ---- 10c. Save model results CSVs ----
write.csv(mixed_results$model_results$performance,
          file.path(output_tab_dir, "model_performance_detailed_LightGBM.csv"), row.names = FALSE)
write.csv(mixed_results$results_summary,
          file.path(output_tab_dir, "counterfactual_summary_LightGBM.csv"), row.names = FALSE)
write.csv(mixed_results$counterfactual_results$results,
          file.path(output_tab_dir, "counterfactual_detailed_LightGBM.csv"), row.names = FALSE)
cat("Model performance CSVs saved to", output_tab_dir, "\n")

#===================== 11. LightGBM feature importance plot =====================#

cat("Generating LightGBM feature importance plot...\n")
lightgbm_model <- mixed_results$model_results$lightgbm_model
if (!is.null(lightgbm_model)) {
  tryCatch({
    importance_matrix <- lgb.importance(lightgbm_model)
    if (!is.null(importance_matrix) && nrow(importance_matrix) > 0) {
      p_importance <- ggplot(importance_matrix,
                             aes(x = reorder(Feature, Gain), y = Gain)) +
        geom_col(fill = "steelblue", alpha = 0.8) +
        coord_flip() +
        labs(
          title = "LightGBM model feature importance",
          subtitle = "Based on information gain (Gain)",
          x = "Feature", y = "Information gain"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          axis.text.y = element_text(size = 10),
          panel.grid.major.y = element_blank()
        )

      ggsave(file.path(output_fig_dir, "LightGBM_feature_importance.png"),
             p_importance, width = 10, height = 8, dpi = 300)
      cat("✓ LightGBM feature importance plot saved successfully\n")
      print(p_importance)
    } else cat("Cannot retrieve LightGBM feature importance\n")
  }, error = function(e) cat("✗ LightGBM feature importance plot failed:", e$message, "\n"))
} else cat("LightGBM model not found, skipping feature importance plot\n")

#===================== 12. LightGBM SHAP analysis & combined Figure 5 =====================#

cat("\n=== Starting LightGBM SHAP analysis (full sample) ===\n")
shap_results <- create_lightgbm_shap_analysis_full(mixed_results$model_results)

cat("Creating SHAP explainability plot (density-weighted beeswarm)...\n")
shap_explainability_plot <- create_shap_explainability_plot(shap_results, mixed_results$model_results, output_fig_dir)

cat("Creating feature importance vs SHAP combined plot...\n")
combined_imp_shap_plot <- create_combined_importance_shap_plot(shap_results, mixed_results$model_results, output_fig_dir)

cat("Creating advanced SHAP summary plot...\n")
advanced_shap_plot <- create_advanced_shap_summary_plot(shap_results, mixed_results$model_results, output_fig_dir)

# ---- 12a. Save SHAP analysis outputs ----
if (!is.null(shap_results)) {
  cat("SHAP analysis succeeded, saving individual figures and data...\n")
  shap_dir <- file.path(output_fig_dir, "SHAP_analysis")
  if (!dir.exists(shap_dir)) dir.create(shap_dir, recursive = TRUE)

  tryCatch({
    ggsave(file.path(shap_dir, "LightGBM_SHAP_summary.png"),
           shap_results$shap_summary, width = 12, height = 10, dpi = 300)
    cat("✓ SHAP summary plot saved\n")
  }, error = function(e) cat("✗ SHAP summary plot save failed:", e$message, "\n"))

  tryCatch({
    ggsave(file.path(shap_dir, "LightGBM_SHAP_waterfall.png"),
           shap_results$waterfall, width = 14, height = 10, dpi = 300)
    cat("✓ SHAP waterfall plot saved\n")
  }, error = function(e) cat("✗ SHAP waterfall plot save failed:", e$message, "\n"))

  tryCatch({
    write.csv(as.data.frame(shap_results$shap_values),
              file.path(shap_dir, "SHAP_values_LightGBM.csv"), row.names = TRUE)
    cat("✓ SHAP values CSV saved\n")
  }, error = function(e) cat("✗ SHAP values CSV save failed:", e$message, "\n"))

  shap_report_content <- paste(
    "=== LightGBM model SHAP analysis report ===\n\n",
    "Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
    "1. Overview\n",
    "SHAP (SHapley Additive exPlanations) is used to explain ML model predictions.\n\n",
    "2. Key indicators\n",
    "Baseline prediction:", round(shap_results$baseline, 4), "\n",
    "Number of samples:", nrow(shap_results$shap_values), "\n",
    "Number of features:", ncol(shap_results$shap_values), "\n\n",
    "3. File list:\n",
    "- LightGBM_SHAP_summary.png\n",
    "- LightGBM_SHAP_waterfall.png\n",
    "- SHAP_values_LightGBM.csv\n",
    "- SHAP_explainability_beeswarm.png\n",
    "- LightGBM_importance_vs_SHAP.png\n",
    "- LightGBM_SHAP_advanced_summary.png\n",
    sep = ""
  )
  writeLines(shap_report_content,
             file.path(shap_dir, "SHAP_analysis_report_LightGBM.txt"))
  cat("✓ SHAP analysis report saved\n")

  print(shap_results$shap_summary)
  print(shap_results$waterfall)
} else {
  cat("SHAP analysis failed, only feature importance plot generated as alternative.\n")
}

# ---- 12b. Assemble Figure 5: counterfactual (top) + SHAP beeswarm (bottom) ----
cat("\n=== Assembling combined Figure 5 ===\n")

if (!is.null(shap_explainability_plot)) {
  # Remove internal title from SHAP plot to avoid duplication with manuscript caption
  fig5_shap <- shap_explainability_plot +
    labs(title = NULL, subtitle = NULL) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )

  # Wrap each part as an atomic element, then stack vertically
  fig5_combined <- (wrap_elements(combined_plot_main) / wrap_elements(fig5_shap)) +
    plot_layout(heights = c(1.2, 0.8))

  ggsave(file.path(output_fig_dir, "Figure5_combined_counterfactual_SHAP.png"),
         fig5_combined, width = 20, height = 18, dpi = 300)
  cat("✓ Figure 5 (combined counterfactual + SHAP) saved\n")

  # Also save the standalone 4-panel counterfactual as a component
  ggsave(file.path(output_fig_dir, "BAA_improvement_combined_LightGBM.png"),
         combined_plot_main, width = 18, height = 12, dpi = 300)
  cat("✓ Standalone 4-panel counterfactual saved\n")
} else {
  # Fallback: save 4-panel alone if SHAP failed
  ggsave(file.path(output_fig_dir, "Figure5_counterfactual_only.png"),
         combined_plot_main, width = 18, height = 12, dpi = 300)
  cat("⚠ SHAP plot unavailable; Figure 5 saved as counterfactual-only\n")
}

#===================== 13. Supplementary figure: by age group × improvement type boxplot =====================#

cat("\n=== Generating by age group × exposure improvement type & Table X (ML) ===\n")

# 1) Prepare data for boxplot and Table X ---------------------------

cf_res <- mixed_results$counterfactual_results$results

# Age group labels (English)
age_breaks  <- c(0, 60, 70, 80, 100)
age_labels  <- c("≤60 years", "61-70 years", "71-80 years", ">80 years")

cf_res2 <- cf_res %>%
  mutate(
    AgeGroup = cut(age,
                   breaks = age_breaks,
                   labels = age_labels,
                   right  = TRUE),
    # Exposure improvement pattern:
    #   Both improved  : both EVI and O3 improved
    #   Only O3 improved : only O3 improved
    ImprovementType = dplyr::case_when(
      evi_intervention & o3_intervention ~ "Both improved",
      !evi_intervention & o3_intervention ~ "Only O3 improved",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(ImprovementType))

cf_res2$ImprovementType <- factor(
  cf_res2$ImprovementType,
  levels = c("Both improved", "Only O3 improved")
)

# Add "All participants"
cf_all <- cf_res2 %>%
  mutate(AgeGroup = factor(
    "All participants",
    levels = c(age_labels, "All participants")
  ))

box_df <- bind_rows(cf_res2, cf_all)

box_df$AgeGroup <- factor(
  box_df$AgeGroup,
  levels = c("≤60 years", "61-70 years", "71-80 years", ">80 years", "All participants")
)

# 2) Table X: mean, SE, and 95%CI by age group × improvement type -------------

table_age <- cf_res2 %>%
  group_by(AgeGroup, ImprovementType) %>%
  summarise(
    N      = n(),
    Coeff  = mean(BAA_Reduction),
    StdErr = sd(BAA_Reduction) / sqrt(N),
    CI_lower = Coeff - 1.96 * StdErr,
    CI_upper = Coeff + 1.96 * StdErr,
    .groups = "drop"
  )

table_all <- cf_all %>%
  group_by(AgeGroup, ImprovementType) %>%
  summarise(
    N      = n(),
    Coeff  = mean(BAA_Reduction),
    StdErr = sd(BAA_Reduction) / sqrt(N),
    CI_lower = Coeff - 1.96 * StdErr,
    CI_upper = Coeff + 1.96 * StdErr,
    .groups = "drop"
  )

table_X <- bind_rows(table_age, table_all) %>%
  mutate(
    AgeGroup = factor(
      AgeGroup,
      levels = c("≤60 years", "61-70 years", "71-80 years", ">80 years", "All participants")
    ),
    ImprovementType = factor(
      ImprovementType,
      levels = c("Both improved", "Only O3 improved")
    )
  ) %>%
  arrange(AgeGroup, ImprovementType)

cat("\n--- Table X: BAA improvement estimates based on ML counterfactual prediction ---\n")
print(table_X)

write.csv(
  table_X,
  file.path(output_tab_dir, "TableX_BAA_improvement_by_age_and_pattern_ML.csv"),
  row.names = FALSE
)
cat("✓ Table X saved as: TableX_BAA_improvement_by_age_and_pattern_ML.csv\n")

# 3) Generate publication-quality boxplot -------------------------------

plot_df <- box_df %>%
  filter(ImprovementType %in% c("Both improved", "Only O3 improved"))

# Prepare data for n labels
n_df <- plot_df %>%
  group_by(AgeGroup, ImprovementType) %>%
  summarise(N = n(), .groups = "drop")

# Max value per age group, used to place n labels
y_range <- range(plot_df$BAA_Reduction, na.rm = TRUE)
y_top   <- y_range[2] + 0.06 * diff(y_range)  # slightly higher globally

n_df <- n_df %>%
  mutate(y_pos = y_top,
         label = paste0("n = ", N))

# Main plot
p_box_age_imp <- ggplot(
  plot_df,
  aes(x = ImprovementType, y = BAA_Reduction)
) +
  # Scatter layer (grey, with slight jitter)
  geom_jitter(
    width = 0.18,
    alpha = 0.30,
    size  = 0.7,
    color = "grey40"
  ) +
  # Boxplot layer, different fill colours per improvement pattern
  geom_boxplot(
    aes(fill = ImprovementType),
    width = 0.5,
    alpha = 0.75,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.4
  ) +
  # n labels
  geom_text(
    data = n_df,
    aes(x = ImprovementType, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 4.2,
    color = "navyblue",
    fontface = "bold"
  ) +
  facet_grid(
    ~ AgeGroup,
    switch = "x"          # Age group strip placed at bottom
  ) +
  coord_cartesian(ylim = c(y_range[1], y_top + 0.04 * diff(y_range))) +
  scale_fill_manual(
    name   = "Improvement pattern",
    values = c("Both improved" = "#F8766D",   # red palette
               "Only O3 improved" = "#00BFC4") # blue-green palette
  ) +
  labs(
    title    = "BAA improvement by age group and exposure improvement pattern",
    y        = "BAA improvement (years)",
    x        = NULL,
    caption  = "Positive values indicate improvement (lower biological age). 'All participants' refers to 45 < age < 90."
  ) +
  theme_classic(base_size = 16) +
  theme(
    # Title and subtitle
    plot.title = element_text(
      family = "serif", face = "bold",
      size = 22, hjust = 0.5, margin = margin(b = 10)
    ),
    plot.caption = element_text(
      size = 11, hjust = 0.5,
      margin = margin(t = 12)
    ),
    # Axes
    axis.title.y = element_text(
      size = 18, face = "bold",
      margin = margin(r = 10)
    ),
    axis.text.y  = element_text(size = 14),
    axis.text.x  = element_text(
      size = 14, angle = 40,
      hjust = 1, vjust = 1
    ),
    # Facet strip at bottom, representing age groups
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(
      family = "serif", face = "bold",
      size = 16, colour = "navy"
    ),
    # Grid lines (horizontal dashed only)
    panel.grid.major.y = element_line(
      colour = "grey90", linetype = "dashed"
    ),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    # Legend
    legend.position = "top",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14)
  )

print(p_box_age_imp)

ggsave(
  file.path(output_fig_dir, "BAA_improvement_boxplot_by_age_ML.png"),
  p_box_age_imp,
  width = 12, height = 7, dpi = 300
)
cat("✓ Boxplot saved as: BAA_improvement_boxplot_by_age_ML.png\n")

# 4) Export boxplot summary statistics (median, quartiles, etc., for appendix) ---------

box_stats <- plot_df %>%
  group_by(AgeGroup, ImprovementType) %>%
  summarise(
    N      = n(),
    Median = median(BAA_Reduction),
    Q1     = quantile(BAA_Reduction, 0.25),
    Q3     = quantile(BAA_Reduction, 0.75),
    Min    = min(BAA_Reduction),
    Max    = max(BAA_Reduction),
    .groups = "drop"
  ) %>%
  arrange(AgeGroup, ImprovementType)

write.csv(
  box_stats,
  file.path(output_tab_dir, "BAA_improvement_boxplot_stats_ML.csv"),
  row.names = FALSE
)
cat("✓ Boxplot stats saved as: BAA_improvement_boxplot_stats_ML.csv\n")

#===================== 14. Final file listing =====================#

cat("\n=== Final file list ===\n")
for (d in c(output_fig_dir, output_tab_dir)) {
  files <- list.files(d, full.names = FALSE, recursive = TRUE)
  cat("Directory", d, "contains", length(files), "files:\n")
  for (f in files) {
    fp <- file.path(d, f)
    sz <- file.info(fp)$size
    cat("  - ", f, " (", format(sz, big.mark = ","), " bytes)\n", sep = "")
  }
}

cat("\n=== All analysis and figure generation completed ===\n")


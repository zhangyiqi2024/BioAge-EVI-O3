############################################################
# Common Section: Load Packages & Read Data
############################################################

# Load required packages
set.seed(20240101)

library(haven)      # Read Stata data
library(ggplot2)    # Plotting
library(dplyr)      # Data cleaning
library(purrr)      # Functional programming
library(mgcv)       # Generalized additive model
library(splines)    # Spline function
library(boot)       # Bootstrap
library(MASS)       # 多元统计分析
library(patchwork)  # Combine multiple ggplot plots

# Set working directory and read data
data <- read_dta("DATA/analysis_ready/analysis_dataset.dta")

# Back up raw data; EVI and O3 processed separately
data_raw <- data

# Check data frame structure and variable names
cat("Raw data frame dimensions:", dim(data_raw), "\n")
cat("变量名:\n")
print(names(data_raw))

cat("First few rows of data:\n")
print(head(data_raw))


############################################################
# Part 1: EVI Exposure -- Method 3: GPSM Spline Regression Model
############################################################

cat("\n==============================\n")
cat("EVI 暴露：GPSM 样条回归模型\n")
cat("==============================\n")

data_evi <- data_raw

# Check if required variables exist（EVI 相关）
required_vars_evi <- c(
  "BAA_PhenoAge", "EVI_1KM_mean4y", "age", "sex", "education", 
  "hh_consumption_total", "retirement", "pension_insurance", "smoking_current", "alcohol_current", 
  "TMP", "PER", "SNOW", "DEM", "lat", 
  "bmi", "hypertension", "dyslipidemia", "year"
)

existing_vars_evi <- required_vars_evi[required_vars_evi %in% names(data_evi)]
missing_vars_evi  <- required_vars_evi[!required_vars_evi %in% names(data_evi)]

cat("EVI分析 - 存在的变量:", paste(existing_vars_evi, collapse = ", "), "\n")
cat("EVI分析 - 缺失的变量:", paste(missing_vars_evi, collapse = ", "), "\n")

# Stop if key variables are missing
if (!"BAA_PhenoAge" %in% names(data_evi)) {
  stop("Dependent variable BAA_PhenoAge not found in data（EVI分析）")
}
if (!"EVI_1KM_mean4y" %in% names(data_evi)) {
  stop("Independent variable EVI_1KM_mean4y not found in data（EVI分析）")
}

# Data preprocessing
if("age" %in% names(data_evi)) {
  data_evi <- data_evi[data_evi$age > 45 & data_evi$age < 90, ]
}

# 创建对数hh_consumption_total
if("hh_consumption_total" %in% names(data_evi)) {
  data_evi$ln_hh_consumption_total <- ifelse(data_evi$hh_consumption_total <= 0, NA, log(data_evi$hh_consumption_total))
  existing_vars_evi <- c(existing_vars_evi, "ln_hh_consumption_total")
}

# Select existing variables
df_evi <- data_evi[, names(data_evi) %in% existing_vars_evi]

# Remove NA
df_evi <- df_evi[complete.cases(df_evi), ]

# Remove rows with infinite values (numeric columns only)
numeric_cols_evi <- sapply(df_evi, is.numeric)
if(any(numeric_cols_evi)) {
  infinite_rows_evi <- apply(df_evi[, numeric_cols_evi, drop = FALSE], 1,
                             function(x) any(is.infinite(x)))
  df_evi <- df_evi[!infinite_rows_evi, ]
}

# Convert haven labels to factors
for(col in names(df_evi)) {
  if(inherits(df_evi[[col]], "labelled")) {
    df_evi[[col]] <- as_factor(df_evi[[col]])
  }
}

cat("EVI分析 - 处理后的数据维度:", dim(df_evi), "\n")
cat("EVI分析 - 处理后的变量名:\n")
print(names(df_evi))

if (nrow(df_evi) == 0) {
  stop("EVI分析 - 数据处理后没有剩余观测值")
}

############################################################
# Climate PCA（EVI分析使用）
############################################################

cat("EVI分析 - Climate PCA...\n")

climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")
climate_data_evi <- df_evi[, climate_vars]

if(any(is.na(climate_data_evi))) {
  cat("Climate variables have missing values; removing affected rows（EVI分析）\n")
  complete_cases_evi <- complete.cases(climate_data_evi)
  climate_data_evi <- climate_data_evi[complete_cases_evi, ]
  df_evi <- df_evi[complete_cases_evi, ]
}

pca_result_evi <- prcomp(climate_data_evi, scale. = TRUE)
cat("EVI分析 - PCA结果摘要:\n")
print(summary(pca_result_evi))

climate_pcs_evi <- as.data.frame(pca_result_evi$x[, 1:5])
colnames(climate_pcs_evi) <- paste0("climate_pc", 1:5)

df_evi <- cbind(df_evi, climate_pcs_evi)

# Convert categorical variables to factors
df_model_evi <- df_evi
factor_vars_evi <- c("sex", "education", "retirement", "pension_insurance", 
                     "smoking_current", "alcohol_current", "hypertension", "dyslipidemia", "year")

for (var in factor_vars_evi) {
  if (var %in% names(df_model_evi)) {
    df_model_evi[[var]] <- as.factor(df_model_evi[[var]])
  }
}

cat("EVI分析 - 最终数据集维度:", dim(df_model_evi), "\n")
cat("EVI分析 - 数据摘要:\n")
summary(df_model_evi)

############################################################
# Method 3: GPSM Spline Regression Model (EVI)
############################################################

cat("\n使用GPSM样条回归模型进行EVI剂量反应分析...\n")

# Control variables
control_vars_evi <- c(
  "age", "sex", "education", "ln_hh_consumption_total", "retirement", 
  "pension_insurance", "smoking_current", "alcohol_current", 
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4", "climate_pc5",
  "bmi", "hypertension", "dyslipidemia", "year"
)

existing_control_vars_evi <- control_vars_evi[control_vars_evi %in% names(df_model_evi)]
cat("EVI分析 - 使用的Control variables:", paste(existing_control_vars_evi, collapse = ", "), "\n")

# EVI Min-Max normalization
df_model_evi$mmx_EVI_1KM_mean4y <- (df_model_evi$EVI_1KM_mean4y -
                                      min(df_model_evi$EVI_1KM_mean4y)) /
  (max(df_model_evi$EVI_1KM_mean4y) - min(df_model_evi$EVI_1KM_mean4y))

# Quantile grouping (5 quintiles)
quantiles_evi <- quantile(df_model_evi$mmx_EVI_1KM_mean4y, probs = seq(0, 1, 0.2))
df_model_evi$cut_better <- as.numeric(cut(df_model_evi$mmx_EVI_1KM_mean4y,
                                          breaks = quantiles_evi,
                                          include.lowest = TRUE))

cat("EVI分析 - 分位数分组效果检查:\n")
group_stats_evi <- by(df_model_evi$mmx_EVI_1KM_mean4y, df_model_evi$cut_better,
                      function(x) c(mean = mean(x), count = length(x), sd = sd(x)))
print(group_stats_evi)

# cut 变量
df_model_evi$cut <- NA
df_model_evi$cut[df_model_evi$cut_better == 1] <- 0.2
df_model_evi$cut[df_model_evi$cut_better == 2] <- 0.4
df_model_evi$cut[df_model_evi$cut_better == 3] <- 0.6
df_model_evi$cut[df_model_evi$cut_better == 4] <- 0.8
df_model_evi$cut[df_model_evi$cut_better == 5] <- 1.0

cat("EVI分析 - 最终GPS分组(分位数法):\n")
final_group_stats_evi <- by(df_model_evi$mmx_EVI_1KM_mean4y, df_model_evi$cut,
                            function(x) c(mean = mean(x), count = length(x)))
print(final_group_stats_evi)

# Estimate GPS via linear regression
cat("EVI分析 - 使用Estimate GPS via linear regression...\n")

if(length(existing_control_vars_evi) > 0) {
  formula_str_evi <- paste("mmx_EVI_1KM_mean4y ~",
                           paste(existing_control_vars_evi, collapse = " + "))
} else {
  formula_str_evi <- "mmx_EVI_1KM_mean4y ~ 1"
}
cat("EVI分析 - 使用的GPS模型公式:", formula_str_evi, "\n")

gps_model_evi <- lm(as.formula(formula_str_evi), data = df_model_evi)

df_model_evi$predicted_EVI <- predict(gps_model_evi)
df_model_evi$residuals <- residuals(gps_model_evi)
residual_sd_evi <- sd(df_model_evi$residuals)

df_model_evi$gps <- dnorm(df_model_evi$mmx_EVI_1KM_mean4y,
                          mean = df_model_evi$predicted_EVI,
                          sd = residual_sd_evi)

df_model_evi$gps_weights <- 1 / df_model_evi$gps

# Trim weights
weight_cutoff_evi <- quantile(df_model_evi$gps_weights, 0.99, na.rm = TRUE)
df_model_evi$gps_weights_trimmed <- pmin(df_model_evi$gps_weights, weight_cutoff_evi)

# Weighted spline regression
cat("EVI分析 - 拟合GPSM样条回归模型...\n")
gpsm_spline_model_evi <- lm(
  BAA_PhenoAge ~ ns(mmx_EVI_1KM_mean4y, df = 2),
  data   = df_model_evi,
  weights = gps_weights_trimmed
)

cat("EVI分析 - GPSM样条模型摘要:\n")
summary(gpsm_spline_model_evi)

# Dose-response curve prediction
x_range_gpsm_evi <- seq(min(df_model_evi$mmx_EVI_1KM_mean4y),
                        max(df_model_evi$mmx_EVI_1KM_mean4y),
                        length.out = 100)

pred_gpsm_data_evi <- data.frame(
  mmx_EVI_1KM_mean4y = x_range_gpsm_evi
)

pred_gpsm_evi <- predict(gpsm_spline_model_evi, newdata = pred_gpsm_data_evi, se.fit = TRUE)

dose_response_gpsm <- data.frame(
  EVI_exposure = x_range_gpsm_evi *
    (max(df_model_evi$EVI_1KM_mean4y) - min(df_model_evi$EVI_1KM_mean4y)) +
    min(df_model_evi$EVI_1KM_mean4y),
  predicted_BAA = pred_gpsm_evi$fit,
  lower = pred_gpsm_evi$fit - 1.96 * pred_gpsm_evi$se.fit,
  upper = pred_gpsm_evi$fit + 1.96 * pred_gpsm_evi$se.fit
)

# Plotting：全部改为英文
p_gpsm <- ggplot(dose_response_gpsm, aes(x = EVI_exposure, y = predicted_BAA)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.3, fill = "orange") +
  geom_line(size = 1.2, color = "darkorange") +
  labs(
    title = "Dose–response curve for EVI",
    subtitle = "Estimated using a generalized propensity score model with spline terms",
    x = "EVI exposure (EVI_1KM_mean4y)",
    y = "Predicted phenotypic age (BAA_PhenoAge)",
    caption = "Shaded area indicates 95% confidence interval; GPS weighting adjusts for confounding"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

print(p_gpsm)

# Save EVI-GPSM plot separately
ggsave("output/figures/剂量反应曲线_GPSM样条_EVI.png",
       p_gpsm, width = 10, height = 6, dpi = 300)


############################################################
# Part 2: O3 Exposure -- Method 1: Simplified GAM Model
############################################################

cat("\n==============================\n")
cat("O3 暴露：简化GAM模型\n")
cat("==============================\n")

data_o3 <- data_raw

# Check if required variables exist（O3 相关）
required_vars_o3 <- c(
  "BAA_PhenoAge", "O3_1KM_mean4y", "age", "sex", "education", 
  "hh_consumption_total", "retirement", "pension_insurance", "smoking_current", "alcohol_current", 
  "hobbies", "TMP", "PER", "SNOW", "DEM", "lat", 
  "bmi", "hypertension", "dyslipidemia", "year"
)

existing_vars_o3 <- required_vars_o3[required_vars_o3 %in% names(data_o3)]
missing_vars_o3  <- required_vars_o3[!required_vars_o3 %in% names(data_o3)]

cat("O3分析 - 存在的变量:", paste(existing_vars_o3, collapse = ", "), "\n")
cat("O3分析 - 缺失的变量:", paste(missing_vars_o3, collapse = ", "), "\n")

# Stop if key variables are missing
if (!"BAA_PhenoAge" %in% names(data_o3)) {
  stop("Dependent variable BAA_PhenoAge not found in data（O3分析）")
}
if (!"O3_1KM_mean4y" %in% names(data_o3)) {
  stop("Independent variable O3_1KM_mean4y not found in data（O3分析）")
}

# Data preprocessing
if("age" %in% names(data_o3)) {
  data_o3 <- data_o3[data_o3$age > 45 & data_o3$age < 90, ]
}

# 创建对数hh_consumption_total
if("hh_consumption_total" %in% names(data_o3)) {
  data_o3$ln_hh_consumption_total <- ifelse(data_o3$hh_consumption_total <= 0, NA, log(data_o3$hh_consumption_total))
  existing_vars_o3 <- c(existing_vars_o3, "ln_hh_consumption_total")
}

# Select existing variables
df_o3 <- data_o3[, names(data_o3) %in% existing_vars_o3]

# Remove NA
df_o3 <- df_o3[complete.cases(df_o3), ]

# Remove rows with infinite values (numeric columns only)
numeric_cols_o3 <- sapply(df_o3, is.numeric)
if(any(numeric_cols_o3)) {
  infinite_rows_o3 <- apply(df_o3[, numeric_cols_o3, drop = FALSE], 1,
                            function(x) any(is.infinite(x)))
  df_o3 <- df_o3[!infinite_rows_o3, ]
}

# Convert haven labels to factors
for(col in names(df_o3)) {
  if(inherits(df_o3[[col]], "labelled")) {
    df_o3[[col]] <- as_factor(df_o3[[col]])
  }
}

cat("O3分析 - 处理后的数据维度:", dim(df_o3), "\n")
cat("O3分析 - 处理后的变量名:\n")
print(names(df_o3))

if (nrow(df_o3) == 0) {
  stop("O3分析 - 数据处理后没有剩余观测值")
}

# Convert categorical variables to factors
df_model_o3 <- df_o3
factor_vars_o3 <- c("sex", "education", "retirement", "pension_insurance", 
                    "smoking_current", "alcohol_current", "hypertension", "dyslipidemia", "year")

for (var in factor_vars_o3) {
  if (var %in% names(df_model_o3)) {
    df_model_o3[[var]] <- as.factor(df_model_o3[[var]])
  }
}

cat("O3分析 - 最终数据集维度:", dim(df_model_o3), "\n")
cat("O3分析 - 数据摘要:\n")
summary(df_model_o3)

############################################################
# Method 1: Simplified GAM Model (O3)
############################################################

cat("使用简化版GAM模型进行O3剂量反应分析...\n")

# GAM smooth terms
gam_terms_o3 <- c("s(O3_1KM_mean4y, k = 50)")

if("age" %in% names(df_model_o3))        gam_terms_o3 <- c(gam_terms_o3, "s(age, k = 5)")
if("bmi" %in% names(df_model_o3)) gam_terms_o3 <- c(gam_terms_o3, "s(bmi, k = 5)")

# Linear terms
linear_vars_o3 <- c(
  "sex", "education", "ln_hh_consumption_total", "retirement", "pension_insurance", 
  "smoking_current", "alcohol_current", "hobbies", "hypertension", "dyslipidemia",
  "TMP", "PER", "SNOW", "DEM", "lat"
)
existing_linear_vars_o3 <- linear_vars_o3[linear_vars_o3 %in% names(df_model_o3)]

# Build GAM formula
gam_formula_o3 <- as.formula(
  paste("BAA_PhenoAge ~", paste(c(gam_terms_o3, existing_linear_vars_o3),
                                collapse = " + "))
)

cat("O3分析 - GAM模型公式:", deparse(gam_formula_o3), "\n")

# Fit GAM model
gam_model_o3 <- gam(
  gam_formula_o3,
  data   = df_model_o3,
  method = "REML"
)

cat("O3分析 - GAM模型摘要:\n")
summary(gam_model_o3)

# Dose-response prediction
x_range_o3 <- seq(min(df_model_o3$O3_1KM_mean4y),
                  max(df_model_o3$O3_1KM_mean4y),
                  length.out = 100)

pred_data_o3 <- df_model_o3[1, ]
pred_data_o3 <- pred_data_o3[rep(1, length(x_range_o3)), ]

pred_data_o3$O3_1KM_mean4y <- x_range_o3

# Continuous covariates set to mean
continuous_vars_o3 <- c("age", "ln_hh_consumption_total", "TMP", "PER", "SNOW", "DEM", "lat", "bmi")
for (var in continuous_vars_o3) {
  if (var %in% names(df_model_o3) && is.numeric(df_model_o3[[var]])) {
    pred_data_o3[[var]] <- mean(df_model_o3[[var]], na.rm = TRUE)
  }
}

# Categorical covariates set to mode
categorical_vars_o3 <- c("sex", "education", "retirement", "pension_insurance", 
                         "smoking_current", "alcohol_current", "hobbies", "hypertension", "dyslipidemia")
for (var in categorical_vars_o3) {
  if (var %in% names(df_model_o3) && is.factor(df_model_o3[[var]])) {
    tab <- table(df_model_o3[[var]])
    pred_data_o3[[var]] <- names(tab)[which.max(tab)]
  }
}

pred_o3 <- predict(gam_model_o3, newdata = pred_data_o3, se.fit = TRUE)

dose_response_o3 <- data.frame(
  O3_exposure  = x_range_o3,
  predicted_BAA = pred_o3$fit,
  lower        = pred_o3$fit - 1.96 * pred_o3$se.fit,
  upper        = pred_o3$fit + 1.96 * pred_o3$se.fit
)

# Plotting：全部改为英文
p <- ggplot(dose_response_o3, aes(x = O3_exposure, y = predicted_BAA)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.3, fill = "steelblue") +
  geom_line(size = 1.2, color = "darkblue") +
  labs(
    title = "Dose–response curve for O3",
    subtitle = "Estimated using a generalized additive model (GAM) with a smooth term for O3",
    x = "O3 exposure (O3_1KM_mean4y)",
    y = "Predicted phenotypic age (BAA_PhenoAge)",
    caption = "Shaded area indicates 95% confidence interval"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

print(p)

# 单独保存 O3-GAM 图（文件名保持中文）
ggsave("output/figures/剂量反应曲线_GAM_O3.png",
       p, width = 10, height = 6, dpi = 300)


############################################################
# Calculate significant thresholds and CIs for EVI and O3
# （Using lowest exposure as reference; find first point where 95% CI excludes 0）
############################################################

# EVI 阈值计算
ref_index_evi <- 1  # 最低暴露对应的索引
ref_mu_evi <- pred_gpsm_evi$fit[ref_index_evi]
ref_se_evi <- pred_gpsm_evi$se.fit[ref_index_evi]

effect_evi <- pred_gpsm_evi$fit - ref_mu_evi
effect_se_evi <- sqrt(pred_gpsm_evi$se.fit^2 + ref_se_evi^2)
effect_lower_evi <- effect_evi - 1.96 * effect_se_evi
effect_upper_evi <- effect_evi + 1.96 * effect_se_evi

# Determine effect direction (sign at highest exposure)
direction_evi <- sign(effect_evi[length(effect_evi)])
if (direction_evi >= 0) {
  sig_idx_evi <- which(effect_lower_evi > 0)
} else {
  sig_idx_evi <- which(effect_upper_evi < 0)
}

if (length(sig_idx_evi) > 0) {
  threshold_index_evi <- sig_idx_evi[1]
  threshold_exposure_evi <- dose_response_gpsm$EVI_exposure[threshold_index_evi]
  threshold_effect_evi   <- effect_evi[threshold_index_evi]
  threshold_lower_evi    <- effect_lower_evi[threshold_index_evi]
  threshold_upper_evi    <- effect_upper_evi[threshold_index_evi]
} else {
  threshold_exposure_evi <- NA
  threshold_effect_evi   <- NA
  threshold_lower_evi    <- NA
  threshold_upper_evi    <- NA
}

# O3 阈值计算
ref_index_o3 <- 1
ref_mu_o3 <- pred_o3$fit[ref_index_o3]
ref_se_o3 <- pred_o3$se.fit[ref_index_o3]

effect_o3 <- pred_o3$fit - ref_mu_o3
effect_se_o3 <- sqrt(pred_o3$se.fit^2 + ref_se_o3^2)
effect_lower_o3 <- effect_o3 - 1.96 * effect_se_o3
effect_upper_o3 <- effect_o3 + 1.96 * effect_se_o3

direction_o3 <- sign(effect_o3[length(effect_o3)])
if (direction_o3 >= 0) {
  sig_idx_o3 <- which(effect_lower_o3 > 0)
} else {
  sig_idx_o3 <- which(effect_upper_o3 < 0)
}

if (length(sig_idx_o3) > 0) {
  threshold_index_o3 <- sig_idx_o3[1]
  threshold_exposure_o3 <- dose_response_o3$O3_exposure[threshold_index_o3]
  threshold_effect_o3   <- effect_o3[threshold_index_o3]
  threshold_lower_o3    <- effect_lower_o3[threshold_index_o3]
  threshold_upper_o3    <- effect_upper_o3[threshold_index_o3]
} else {
  threshold_exposure_o3 <- NA
  threshold_effect_o3   <- NA
  threshold_lower_o3    <- NA
  threshold_upper_o3    <- NA
}

# Organize threshold results and save as .csv
threshold_results <- data.frame(
  exposure              = c("EVI_1KM_mean4y", "O3_1KM_mean4y"),
  threshold_exposure    = c(threshold_exposure_evi, threshold_exposure_o3),
  effect_at_threshold   = c(threshold_effect_evi,   threshold_effect_o3),
  lower_95CI            = c(threshold_lower_evi,    threshold_lower_o3),
  upper_95CI            = c(threshold_upper_evi,    threshold_upper_o3),
  reference_exposure    = c(min(dose_response_gpsm$EVI_exposure),
                            min(dose_response_o3$O3_exposure))
)

write.csv(threshold_results, "output/tables/threshold_results_EVI_O3.csv", row.names = FALSE)


############################################################
# Final step: combine two plots on one canvas (English)
############################################################

combined_plot <- p_gpsm + p +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Dose–response curves of EVI and O3",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold")
    )
  )

print(combined_plot)

ggsave("output/figures/剂量反应曲线_EVI_GPSM_与_O3_GAM_组合图.png",
       combined_plot, width = 12, height = 6, dpi = 300)

cat("\nAll analyses complete!\n",
    "Generated:\n",
    "  1) 剂量反应曲线_GPSM样条_EVI.png\n",
    "  2) 剂量反应曲线_GAM_O3.png\n",
    "  3) 剂量反应曲线_EVI_GPSM_与_O3_GAM_组合图.png\n",
    "  4) 阈值结果_EVI_O3.csv\n")

# =============================================================================
# R code: EVI & O3 single-exposure per-IQR associations with BAA (LMM coefficient estimation)
# Includes dual-exposure models, with NA/NaN/Inf error handling
# =============================================================================

set.seed(20240101)

# Load required R packages (includes scales for number formatting)
required_packages <- c("lme4", "lmerTest", "dplyr", "tidyr",
                       "kableExtra", "haven", "ggplot2", "scales")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# =============================================================================
# 1. Data loading and preprocessing (with enhanced error checking)
# =============================================================================

# Attempt to load data from relative paths
data_loaded <- FALSE
data <- NULL
csv_path <- "DATA/analysis_ready/analysis_dataset.csv"
dta_path <- "DATA/analysis_ready/analysis_dataset.dta"

# Try loading CSV format
if (file.exists(csv_path)) {
  data <- read.csv(csv_path, stringsAsFactors = FALSE)
  data_loaded <- TRUE
  cat("Successfully loaded CSV format data\n")
}

# Try loading Stata format
if (!data_loaded && file.exists(dta_path)) {
  data <- as.data.frame(read_dta(dta_path))
  data_loaded <- TRUE
  cat("Successfully loaded Stata format data\n")
}

# If neither exists, report error
if (!data_loaded) {
  stop("Unable to find data file. Checked paths: ", csv_path, " and ", dta_path)
}

cat("Original data dimensions: ", dim(data), "\n")

# =============================================================================
# 2. Enhanced data cleaning: handle all types of invalid values
# =============================================================================

# Define the list of variables to check
analysis_vars <- c("BAA_PhenoAge", "EVI_1KM_mean4y", "O3_1KM_mean4y",
                   "age", "sex", "id_num",
                   "education", "hh_consumption_total", "retirement", "pension_insurance",
                   "smoking_current", "alcohol_current",
                   "TMP", "PER", "SNOW", "DEM", "lat",
                   "bmi", "hypertension", "dyslipidemia", "year")

# Retain only the variables needed for analysis
analysis_data <- data[intersect(analysis_vars, names(data))]

# Check invalid values for each variable
check_invalid_values <- function(df) {
  invalid_summary <- data.frame(
    Variable = character(),
    NA_count = integer(),
    NaN_count = integer(),
    Inf_count = integer(),
    stringsAsFactors = FALSE
  )

  for (col in names(df)) {
    na_count  <- sum(is.na(df[[col]]))
    nan_count <- sum(is.nan(df[[col]]))
    inf_count <- sum(is.infinite(df[[col]]))

    invalid_summary <- rbind(invalid_summary, data.frame(
      Variable = col,
      NA_count = na_count,
      NaN_count = nan_count,
      Inf_count = inf_count
    ))
  }

  return(invalid_summary)
}

invalid_values <- check_invalid_values(analysis_data)
print("Invalid values per variable:")
print(invalid_values)

# Create log-transformed household consumption variable (with enhanced error handling)
if ("hh_consumption_total" %in% names(analysis_data)) {
  if (any(analysis_data$hh_consumption_total <= 0, na.rm = TRUE)) {
    cat("Found", sum(analysis_data$hh_consumption_total <= 0, na.rm = TRUE),
        "non-positive values, will handle before log transformation\n")
    analysis_data$hh_consumption_total[analysis_data$hh_consumption_total <= 0] <- NA
  }
  analysis_data$ln_hh_consumption_total <- log(analysis_data$hh_consumption_total)
}

# Filter age range (45-90 years)
analysis_data <- analysis_data %>%
  filter(age > 45 & age < 90)

# Convert categorical variables to factors
categorical_vars <- c("sex", "education", "retirement", "pension_insurance",
                      "smoking_current", "alcohol_current", "hypertension", "dyslipidemia", "year")

for (var in categorical_vars) {
  if (var %in% names(analysis_data)) {
    analysis_data[[var]] <- as.factor(analysis_data[[var]])
  }
}

# Remove rows containing any invalid values
analysis_data_clean <- analysis_data
rows_before <- nrow(analysis_data_clean)

for (col in names(analysis_data_clean)) {
  invalid_idx <- is.na(analysis_data_clean[[col]]) |
    is.nan(analysis_data_clean[[col]]) |
    is.infinite(analysis_data_clean[[col]])

  if (sum(invalid_idx) > 0) {
    cat("Variable", col, "has", sum(invalid_idx), "invalid values, will be removed\n")
    analysis_data_clean <- analysis_data_clean[!invalid_idx, ]
  }
}

rows_after <- nrow(analysis_data_clean)
cat("Rows before cleaning:", rows_before, "\n")
cat("Rows after cleaning:", rows_after, "\n")
cat("Rows removed:", rows_before - rows_after, "\n")

final_check <- check_invalid_values(analysis_data_clean)
print("Invalid values per variable after cleaning:")
print(final_check)

if (any(final_check$NA_count > 0 |
        final_check$NaN_count > 0 |
        final_check$Inf_count > 0)) {
  stop("Data still contains invalid values, please check data quality")
} else {
  cat("Data cleaning complete, no invalid values\n")
}

data <- analysis_data_clean

# =============================================================================
# 3. Principal component analysis of climate variables
# =============================================================================

cat("=== Climate variable PCA ===\n")

climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")
climate_data <- data[climate_vars]

if (any(is.na(climate_data))) {
  cat("Climate variables contain missing values, removing incomplete cases\n")
  complete_cases <- complete.cases(climate_data)
  climate_data <- climate_data[complete_cases, ]
  data <- data[complete_cases, ]
}

pca_result <- prcomp(climate_data, scale. = TRUE, center = TRUE)
print(summary(pca_result))

climate_pcs <- as.data.frame(pca_result$x[, 1:5])
colnames(climate_pcs) <- paste0("climate_pc", 1:5)
data <- cbind(data, climate_pcs)

cat("Climate PCA complete, 5 principal components added to data\n")

# =============================================================================
# 4. Calculate IQR (interquartile range) of exposure variables
# =============================================================================

calculate_iqr <- function(x) {
  quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE)
}

iqr_evi <- calculate_iqr(data$EVI_1KM_mean4y)
iqr_o3  <- calculate_iqr(data$O3_1KM_mean4y)

cat("EVI IQR value:", round(iqr_evi, 3), "\n")
cat("O3 IQR value:",  round(iqr_o3, 3),  "\n")

# =============================================================================
# 5. Linear mixed-effects model function (with error handling)
# =============================================================================

run_lmm_model_safe <- function(exposure_var, adjustment_level, data) {
  tryCatch({
    if (adjustment_level == "Unadjusted") {
      formula_str <- paste("BAA_PhenoAge ~", exposure_var, "+ (1 | id_num)")
    } else if (adjustment_level == "Age and sex adjusted") {
      formula_str <- paste("BAA_PhenoAge ~", exposure_var, "+ age + sex + (1 | id_num)")
    } else if (adjustment_level == "Multivariable adjusted") {
      formula_str <- paste("BAA_PhenoAge ~", exposure_var,
                           "+ age + sex + education + ln_hh_consumption_total + retirement + pension_insurance +",
                           "smoking_current + alcohol_current + climate_pc1 + climate_pc2 + climate_pc3 + climate_pc4 + climate_pc5 +",
                           "bmi + hypertension + dyslipidemia + year + (1 | id_num)")
    } else if (adjustment_level == "Dual exposure adjusted") {
      formula_str <- paste("BAA_PhenoAge ~ EVI_1KM_mean4y + O3_1KM_mean4y",
                           "+ age + sex + education + ln_hh_consumption_total + retirement + pension_insurance +",
                           "smoking_current + alcohol_current + climate_pc1 + climate_pc2 + climate_pc3 + climate_pc4 + climate_pc5 +",
                           "bmi + hypertension + dyslipidemia + year + (1 | id_num)")
    }

    vars_in_formula <- all.vars(as.formula(formula_str))
    missing_vars <- setdiff(vars_in_formula, names(data))
    if (length(missing_vars) > 0) {
      warning("Variables missing from formula: ", paste(missing_vars, collapse = ", "))
      return(NULL)
    }

    model <- lmer(as.formula(formula_str), data = data, REML = TRUE)
    model_summary <- summary(model)
    coefficients  <- model_summary$coefficients

    if (adjustment_level == "Dual exposure adjusted") {
      results <- list()

      # EVI
      evi_coef_row <- which(rownames(coefficients) == "EVI_1KM_mean4y")
      if (length(evi_coef_row) > 0) {
        coef_value_evi <- coefficients[evi_coef_row, "Estimate"]
        se_value_evi   <- coefficients[evi_coef_row, "Std. Error"]
        p_value_evi    <- coefficients[evi_coef_row, "Pr(>|t|)"]

        ci_lower_evi <- coef_value_evi - 1.96 * se_value_evi
        ci_upper_evi <- coef_value_evi + 1.96 * se_value_evi

        coef_iqr_evi     <- coef_value_evi * iqr_evi
        ci_lower_iqr_evi <- ci_lower_evi * iqr_evi
        ci_upper_iqr_evi <- ci_upper_evi * iqr_evi

        significance_evi <- ""
        if (p_value_evi < 0.001) {
          significance_evi <- "***"
        } else if (p_value_evi < 0.01) {
          significance_evi <- "**"
        } else if (p_value_evi < 0.05) {
          significance_evi <- "*"
        }

        results$evi <- list(
          exposure        = "EVI_1KM_mean4y",
          adjustment      = adjustment_level,
          coefficient     = coef_value_evi,
          coefficient_iqr = coef_iqr_evi,
          ci_lower        = ci_lower_evi,
          ci_upper        = ci_upper_evi,
          ci_lower_iqr    = ci_lower_iqr_evi,
          ci_upper_iqr    = ci_upper_iqr_evi,
          p_value         = p_value_evi,
          significance    = significance_evi,
          model_converged = TRUE
        )
      }

      # O3
      o3_coef_row <- which(rownames(coefficients) == "O3_1KM_mean4y")
      if (length(o3_coef_row) > 0) {
        coef_value_o3 <- coefficients[o3_coef_row, "Estimate"]
        se_value_o3   <- coefficients[o3_coef_row, "Std. Error"]
        p_value_o3    <- coefficients[o3_coef_row, "Pr(>|t|)"]

        ci_lower_o3 <- coef_value_o3 - 1.96 * se_value_o3
        ci_upper_o3 <- coef_value_o3 + 1.96 * se_value_o3

        coef_iqr_o3     <- coef_value_o3 * iqr_o3
        ci_lower_iqr_o3 <- ci_lower_o3 * iqr_o3
        ci_upper_iqr_o3 <- ci_upper_o3 * iqr_o3

        significance_o3 <- ""
        if (p_value_o3 < 0.001) {
          significance_o3 <- "***"
        } else if (p_value_o3 < 0.01) {
          significance_o3 <- "**"
        } else if (p_value_o3 < 0.05) {
          significance_o3 <- "*"
        }

        results$o3 <- list(
          exposure        = "O3_1KM_mean4y",
          adjustment      = adjustment_level,
          coefficient     = coef_value_o3,
          coefficient_iqr = coef_iqr_o3,
          ci_lower        = ci_lower_o3,
          ci_upper        = ci_upper_o3,
          ci_lower_iqr    = ci_lower_iqr_o3,
          ci_upper_iqr    = ci_upper_iqr_o3,
          p_value         = p_value_o3,
          significance    = significance_o3,
          model_converged = TRUE
        )
      }

      return(results)

    } else {
      exposure_coef_row <- which(rownames(coefficients) == exposure_var)
      if (length(exposure_coef_row) == 0) {
        warning("Variable not found in model: ", exposure_var)
        return(NULL)
      }

      coef_value <- coefficients[exposure_coef_row, "Estimate"]
      se_value   <- coefficients[exposure_coef_row, "Std. Error"]
      p_value    <- coefficients[exposure_coef_row, "Pr(>|t|)"]

      ci_lower <- coef_value - 1.96 * se_value
      ci_upper <- coef_value + 1.96 * se_value

      if (exposure_var == "EVI_1KM_mean4y") {
        coef_iqr      <- coef_value * iqr_evi
        ci_lower_iqr  <- ci_lower * iqr_evi
        ci_upper_iqr  <- ci_upper * iqr_evi
      } else {
        coef_iqr      <- coef_value * iqr_o3
        ci_lower_iqr  <- ci_lower * iqr_o3
        ci_upper_iqr  <- ci_upper * iqr_o3
      }

      significance <- ""
      if (p_value < 0.001) {
        significance <- "***"
      } else if (p_value < 0.01) {
        significance <- "**"
      } else if (p_value < 0.05) {
        significance <- "*"
      }

      return(list(
        exposure        = exposure_var,
        adjustment      = adjustment_level,
        coefficient     = coef_value,
        coefficient_iqr = coef_iqr,
        ci_lower        = ci_lower,
        ci_upper        = ci_upper,
        ci_lower_iqr    = ci_lower_iqr,
        ci_upper_iqr    = ci_upper_iqr,
        p_value         = p_value,
        significance    = significance,
        model_converged = TRUE
      ))
    }
  }, error = function(e) {
    warning("Model fitting failed: ", e$message)
    return(list(
      exposure        = exposure_var,
      adjustment      = adjustment_level,
      coefficient     = NA,
      coefficient_iqr = NA,
      ci_lower        = NA,
      ci_upper        = NA,
      ci_lower_iqr    = NA,
      ci_upper_iqr    = NA,
      p_value         = NA,
      significance    = "",
      model_converged = FALSE,
      error_message   = e$message
    ))
  })
}

# =============================================================================
# 6. Run all models and collect results (including dual-exposure model)
# =============================================================================

exposure_vars <- c("EVI_1KM_mean4y", "O3_1KM_mean4y")
adjustment_levels <- c("Unadjusted", "Age and sex adjusted", "Multivariable adjusted")

results <- list()
model_count <- 0
successful_models <- 0

for (exp_var in exposure_vars) {
  for (adj_level in adjustment_levels) {
    model_count <- model_count + 1
    cat("Running single-exposure model", model_count, ": ", exp_var, " - ", adj_level, "\n")

    model_result <- run_lmm_model_safe(exp_var, adj_level, data)

    if (!is.null(model_result) && model_result$model_converged) {
      results[[paste(exp_var, adj_level, sep = "_")]] <- model_result
      successful_models <- successful_models + 1
    }
  }
}

cat("Running dual-exposure model: EVI and O3 simultaneously adjusted\n")
dual_exposure_result <- run_lmm_model_safe("dual", "Dual exposure adjusted", data)

if (!is.null(dual_exposure_result)) {
  if (!is.null(dual_exposure_result$evi) && dual_exposure_result$evi$model_converged) {
    results[["EVI_Dual exposure adjusted"]] <- dual_exposure_result$evi
    successful_models <- successful_models + 1
  }
  if (!is.null(dual_exposure_result$o3) && dual_exposure_result$o3$model_converged) {
    results[["O3_Dual exposure adjusted"]] <- dual_exposure_result$o3
    successful_models <- successful_models + 1
  }
}

cat("Number of successfully fitted models:", successful_models, "/",
    (length(exposure_vars) * length(adjustment_levels) + 2), "\n")

# =============================================================================
# 7. Compile results and create intermediate table data
# =============================================================================

if (successful_models == 0) {
  stop("All model fitting failed, please check data quality and model specification")
}

results_list <- list()
for (result_name in names(results)) {
  result <- results[[result_name]]
  if (is.list(result) && !is.null(result$model_converged) && result$model_converged) {
    results_list[[result_name]] <- as.data.frame(result)
  }
}

if (length(results_list) == 0) {
  stop("No successfully fitted model results")
}

results_df <- do.call(rbind, results_list)

# Display exposure factor as EVI / O3
results_df$Exposure <- ifelse(grepl("EVI", results_df$exposure), "EVI", "O₃")
results_df$ModelType <- ifelse(results_df$adjustment == "Dual exposure adjusted",
                          "Dual-exposure model", "Single-exposure model")
results_df$DisplayFactor <- ifelse(results_df$adjustment == "Dual exposure adjusted",
                          paste0(results_df$Exposure, " (dual-exposure)"),
                          results_df$Exposure)

# Calculate SE on IQR scale (derived from CI)
results_df <- results_df %>%
  mutate(
    SE_iqr = (ci_upper_iqr - ci_lower_iqr) / (2 * 1.96),
    Exposure_label = ifelse(grepl("EVI", exposure), "EVI_iqr", "O3_iqr"),
    ModelType_en   = ifelse(ModelType == "Dual-exposure model",
                            "Dual-exposure", "Single-exposure"),
    Buffer         = "1KM",
    Adjustment_en  = dplyr::recode(
      adjustment,
      "Unadjusted"            = "Unadjusted",
      "Age and sex adjusted"  = "Age and sex adjusted",
      "Multivariable adjusted" = "Multivariable adjusted",
      "Dual exposure adjusted" = "Dual exposure adjusted"
    )
  )

# Generate display column for main HTML table (English, per-IQR effect + 95% CI)
results_df <- results_df %>%
  mutate(
    Display = ifelse(
      !is.na(coefficient_iqr),
      sprintf("%.3f (%.3f, %.3f)%s",
              coefficient_iqr,
              ci_lower_iqr,
              ci_upper_iqr,
              significance),
      "Model failed"
    )
  )

# =============================================================================
# 8. Create main results HTML table (English headers)
# =============================================================================

wide_results <- results_df %>%
  select(Exposure = Exposure,
         Adjustment = Adjustment_en,
         Display) %>%
  tidyr::pivot_wider(names_from = Adjustment, values_from = Display)

desired_columns <- c("Exposure",
                     "Unadjusted",
                     "Age and sex adjusted",
                     "Multivariable adjusted",
                     "Dual exposure adjusted")
existing_columns <- intersect(desired_columns, names(wide_results))
wide_results <- wide_results[, existing_columns, drop = FALSE]

if (nrow(wide_results) > 0) {
  main_table <- kable(
    wide_results,
    col.names = c("Exposure",
                  "Unadjusted",
                  "Age and sex adjusted",
                  "Multivariable adjusted",
                  "Dual exposure adjusted"),
    align = c("l", "c", "c", "c", "c"),
    caption = "Associations of EVI and O₃ with phenotypic age acceleration per IQR increase (linear mixed-effects models, 95% CI)",
    format = "html"
  ) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = FALSE) %>%
    add_header_above(c(" " = 1, "Adjustment set" = 4)) %>%
    row_spec(0, bold = TRUE) %>%
    column_spec(1, bold = TRUE)

  print(main_table)
  save_kable(main_table, file = "output/tables/iqr_models_results.html")
} else {
  warning("No successfully fitted model results to display")
}

# =============================================================================
# 9. Create visualizations (English titles, 3 decimal places)
# =============================================================================

if (nrow(results_df) > 0) {

  plot_data <- results_df %>%
    mutate(
      Adjustment = factor(
        Adjustment_en,
        levels = c("Unadjusted",
                   "Age and sex adjusted",
                   "Multivariable adjusted",
                   "Dual exposure adjusted")
      ),
      Exposure_en = Exposure  # EVI / O3
    )

  # Figure 1: Per-IQR effect under each adjustment level (point + error bar)
  coefficient_plot <- ggplot(
    plot_data,
    aes(x = Adjustment, y = coefficient_iqr,
        color = DisplayFactor, group = DisplayFactor)
  ) +
    geom_point(position = position_dodge(width = 0.7), size = 3) +
    geom_errorbar(aes(ymin = ci_lower_iqr, ymax = ci_upper_iqr),
                  width = 0.2,
                  position = position_dodge(width = 0.7)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
    labs(
      title    = "Associations of EVI and O₃ with phenotypic age acceleration per IQR increase",
      subtitle = "Linear mixed-effects models under different adjustment sets",
      x        = "Adjustment set",
      y        = "Coefficient per IQR increase in exposure",
      color    = "Exposure term"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      axis.text.x     = element_text(angle = 45, hjust = 1)
    )

  print(coefficient_plot)
  ggsave("output/figures/iqr_models_coefficient_plot.png",
         coefficient_plot,
         width = 5, height = 6, dpi = 300)

  # Figure 2: Group comparison by exposure factor across adjustment sets (bar chart)
  group_plot <- ggplot(
    plot_data,
    aes(x = Exposure_en, y = coefficient_iqr, fill = Adjustment)
  ) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = ci_lower_iqr, ymax = ci_upper_iqr),
                  position = position_dodge(width = 0.8),
                  width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
    labs(
      title    = "Comparative effects of EVI and O₃ on phenotypic age acceleration",
      subtitle = "Per IQR increase in exposure across different adjustment sets",
      x        = "Exposure",
      y        = "Coefficient per IQR increase in exposure",
      fill     = "Adjustment set"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom"
    )

  print(group_plot)
  ggsave("output/figures/iqr_models_group_comparison.png",
         group_plot,
         width = 7, height = 6, dpi = 300)
}

# =============================================================================
# 10. Save detailed results summary (CSV: English column names + per-IQR effects)
# =============================================================================

results_for_csv <- results_df %>%
  mutate(
    Exposure     = Exposure_label,      # EVI_iqr / O3_iqr
    `Model Type` = ModelType_en,        # Single-exposure / Dual-exposure
    Buffer       = Buffer,              # "1KM"
    Adjustment_csv = dplyr::case_when(
      adjustment == "Unadjusted"            ~ "Unadjusted",
      adjustment == "Age and sex adjusted"  ~ "Age and sex",
      adjustment %in% c("Multivariable adjusted", "Dual exposure adjusted") ~ "Multivariable",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(
    factor(Exposure,      levels = c("EVI_iqr", "O3_iqr")),
    factor(`Model Type`,  levels = c("Single-exposure", "Dual-exposure")),
    factor(Adjustment_csv, levels = c("Unadjusted", "Age and sex", "Multivariable"))
  ) %>%
  transmute(
    Exposure,
    `Model Type`,
    Buffer,
    Adjustment      = Adjustment_csv,
    `Coefficient (SE)` = sprintf("%.3f (%.3f)", coefficient_iqr, SE_iqr),
    p_value            = sprintf("%.3f", p_value),
    `95% CI`           = sprintf("[%.3f, %.3f]", ci_lower_iqr, ci_upper_iqr)
  )

write.csv(results_for_csv,
          "output/tables/iqr_models_detailed_results.csv",
          row.names = FALSE)

# =============================================================================
# 11. Output summary information
# =============================================================================

cat("\n")
cat("==================================================\n")
cat("Analysis complete! (including dual-exposure models)\n")
cat("==================================================\n")
cat("Data information:\n")
cat("- Cleaned data dimensions:", nrow(data), "rows,", ncol(data), "columns\n")
cat("- EVI IQR value:", round(iqr_evi, 3), "\n")
cat("- O3 IQR value:",  round(iqr_o3, 3), "\n")
cat("\n")
cat("Model fitting results:\n")
cat("- Single-exposure models:", length(exposure_vars) * length(adjustment_levels), "\n")
cat("- Dual-exposure models: 1 (containing both EVI and O3 coefficients)\n")
cat("- Successfully fitted coefficients:", successful_models, "\n")
cat("\n")
cat("Results files saved:\n")
if (exists("wide_results") && nrow(wide_results) > 0) {
  cat("- HTML table: ../output/tables/iqr_models_results.html\n")
}
if (exists("results_df") && nrow(results_df) > 0) {
  cat("- Combined figure: ../output/figures/iqr_models_coefficient_plot.png\n")
  cat("- Group comparison figure: ../output/figures/iqr_models_group_comparison.png\n")
  cat("- Summary CSV: ../output/tables/iqr_models_detailed_results.csv\n")
}
cat("\n")
cat("Significance legend:\n")
cat("* p < 0.05, ** p < 0.01, *** p < 0.001\n")
cat("==================================================\n")

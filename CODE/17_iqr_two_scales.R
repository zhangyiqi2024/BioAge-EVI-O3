# =============================================================================
# R code: 1KM & 500m buffer, single- and dual-exposure LMM IQR analysis
# =============================================================================

set.seed(20240101)

# Load required R packages
required_packages <- c("lme4", "lmerTest", "dplyr", "tidyr",
                       "kableExtra", "haven", "ggplot2")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# =============================================================================
# 1. Data loading and preprocessing (enhanced error checking)
# =============================================================================

# Try to load data (paths relative to the project root; run from 04_code/)
data_loaded <- FALSE
data <- NULL

# Try loading CSV format
csv_path <- "DATA/analysis_ready/analysis_dataset.csv"
if (file.exists(csv_path)) {
  data <- read.csv(csv_path, stringsAsFactors = FALSE)
  data_loaded <- TRUE
  cat("Successfully loaded CSV format data\n")
}

# Try loading Stata format
if (!data_loaded) {
  dta_path <- "DATA/analysis_ready/analysis_dataset.dta"
  if (file.exists(dta_path)) {
    data <- as.data.frame(read_dta(dta_path))
    data_loaded <- TRUE
    cat("Successfully loaded Stata format data\n")
  }
}

if (!data_loaded) {
  stop("Could not find DATA/analysis_ready/analysis_dataset.dta. ",
       "Run from the project root (the 04_code/ folder).")
}

cat("Original data dimensions:", dim(data), "\n")

# =============================================================================
# 2. Enhanced data cleaning: handle all types of invalid values
# =============================================================================

# Add 500m buffer exposure variables
analysis_vars <- c(
  "BAA_PhenoAge",
  "EVI_1KM_mean4y", "O3_1KM_mean4y",
  "EVI_500m_mean4y", "O3_500m_mean4y",
  "age", "sex", "id_num",
  "education", "hh_consumption_total", "retirement", "pension_insurance",
  "smoking_current", "alcohol_current",
  "TMP", "PER", "SNOW", "DEM", "lat",
  "bmi", "hypertension", "dyslipidemia", "year"
)

analysis_data <- data[intersect(analysis_vars, names(data))]

check_invalid_values <- function(df) {
  invalid_summary <- data.frame(
    Variable = character(),
    Missing_count = integer(),
    NaN_count = integer(),
    Inf_count = integer(),
    stringsAsFactors = FALSE
  )
  for (col in names(df)) {
    na_count  <- sum(is.na(df[[col]]))
    nan_count <- sum(is.nan(df[[col]]))
    inf_count <- sum(is.infinite(df[[col]]))
    invalid_summary <- rbind(
      invalid_summary,
      data.frame(
        Variable = col,
        Missing_count = na_count,
        NaN_count = nan_count,
        Inf_count = inf_count
      )
    )
  }
  invalid_summary
}

invalid_values <- check_invalid_values(analysis_data)
print("Invalid values by variable:")
print(invalid_values)

# Log-transform household consumption
if ("hh_consumption_total" %in% names(analysis_data)) {
  if (any(analysis_data$hh_consumption_total <= 0, na.rm = TRUE)) {
    cat("Found", sum(analysis_data$hh_consumption_total <= 0, na.rm = TRUE),
        "non-positive values, will handle before log transformation\n")
    analysis_data$hh_consumption_total[analysis_data$hh_consumption_total <= 0] <- NA
  }
  analysis_data$ln_hh_consumption_total <- log(analysis_data$hh_consumption_total)
}

# Age filtering
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

# Remove all rows with invalid values
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
print("Invalid values by variable after cleaning:")
print(final_check)

if (any(final_check$Missing_count > 0 |
        final_check$NaN_count   > 0 |
        final_check$Inf_count   > 0)) {
  stop("Data still contains invalid values; please check data quality")
} else {
  cat("Data cleaning complete, no invalid values\n")
}

data <- analysis_data_clean

# =============================================================================
# 3. Climate variable principal component analysis
# =============================================================================

cat("=== Climate variable principal component analysis ===\n")

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

cat("Climate PCA complete, added 5 principal components to data\n")

# =============================================================================
# 4. Calculate IQR for exposure variables: 1KM and 500m
# =============================================================================

calculate_iqr <- function(x) {
  quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE)
}

iqr_evi_1km  <- calculate_iqr(data$EVI_1KM_mean4y)
iqr_o3_1km   <- calculate_iqr(data$O3_1KM_mean4y)
iqr_evi_500  <- calculate_iqr(data$EVI_500m_mean4y)
iqr_o3_500   <- calculate_iqr(data$O3_500m_mean4y)

cat("1KM EVI IQR:",  round(iqr_evi_1km, 3), "\n")
cat("1KM O3  IQR:",  round(iqr_o3_1km, 3), "\n")
cat("500m EVI IQR:", round(iqr_evi_500, 3), "\n")
cat("500m O3  IQR:", round(iqr_o3_500, 3), "\n")

# =============================================================================
# 5. Linear mixed model function (supports different buffers)
# =============================================================================

run_lmm_model_safe <- function(exposure_var, adjustment_level, data,
                               evi_name, o3_name,
                               iqr_evi, iqr_o3,
                               buffer_label) {
  tryCatch({
    # Build formula
    if (adjustment_level == "Unadjusted") {
      formula_str <- paste("BAA_PhenoAge ~", exposure_var, "+ (1 | id_num)")
    } else if (adjustment_level == "Age and sex adjusted") {
      formula_str <- paste("BAA_PhenoAge ~", exposure_var,
                           "+ age + sex + (1 | id_num)")
    } else if (adjustment_level == "Multivariable adjusted") {
      formula_str <- paste("BAA_PhenoAge ~", exposure_var,
                           "+ age + sex + education + ln_hh_consumption_total + retirement + pension_insurance +",
                           "smoking_current + alcohol_current + climate_pc1 + climate_pc2 + climate_pc3 + climate_pc4 + climate_pc5 +",
                           "bmi + hypertension + dyslipidemia + year + (1 | id_num)")
    } else if (adjustment_level == "Dual exposure adjusted") {
      # Dual-exposure model: simultaneously include EVI and O3 (by buffer variable names)
      formula_str <- paste("BAA_PhenoAge ~", evi_name, "+", o3_name,
                           "+ age + sex + education + ln_hh_consumption_total + retirement + pension_insurance +",
                           "smoking_current + alcohol_current + climate_pc1 + climate_pc2 + climate_pc3 + climate_pc4 + climate_pc5 +",
                           "bmi + hypertension + dyslipidemia + year + (1 | id_num)")
    }

    # Check if variables exist
    vars_in_formula <- all.vars(as.formula(formula_str))
    missing_vars <- setdiff(vars_in_formula, names(data))
    if (length(missing_vars) > 0) {
      warning("Missing variables in formula: ", paste(missing_vars, collapse = ", "))
      return(NULL)
    }

    # Fit model
    model <- lmer(as.formula(formula_str), data = data, REML = TRUE)
    model_summary <- summary(model)
    coefficients <- model_summary$coefficients

    # Dual-exposure model
    if (adjustment_level == "Dual exposure adjusted") {
      results <- list()

      # EVI
      evi_coef_row <- which(rownames(coefficients) == evi_name)
      if (length(evi_coef_row) > 0) {
        coef_value_evi <- coefficients[evi_coef_row, "Estimate"]
        se_value_evi   <- coefficients[evi_coef_row, "Std. Error"]
        p_value_evi    <- coefficients[evi_coef_row, "Pr(>|t|)"]

        ci_lower_evi <- coef_value_evi - 1.96 * se_value_evi
        ci_upper_evi <- coef_value_evi + 1.96 * se_value_evi

        coef_iqr_evi      <- coef_value_evi * iqr_evi
        ci_lower_iqr_evi  <- ci_lower_evi * iqr_evi
        ci_upper_iqr_evi  <- ci_upper_evi * iqr_evi

        significance_evi <- ""
        if (p_value_evi < 0.001) {
          significance_evi <- "***"
        } else if (p_value_evi < 0.01) {
          significance_evi <- "**"
        } else if (p_value_evi < 0.05) {
          significance_evi <- "*"
        }

        results$evi <- list(
          exposure        = evi_name,
          adjustment      = adjustment_level,
          coefficient     = coef_value_evi,
          coefficient_iqr = coef_iqr_evi,
          ci_lower        = ci_lower_evi,
          ci_upper        = ci_upper_evi,
          ci_lower_iqr    = ci_lower_iqr_evi,
          ci_upper_iqr    = ci_upper_iqr_evi,
          p_value         = p_value_evi,
          significance    = significance_evi,
          buffer          = buffer_label,
          model_converged = TRUE
        )
      }

      # O3
      o3_coef_row <- which(rownames(coefficients) == o3_name)
      if (length(o3_coef_row) > 0) {
        coef_value_o3 <- coefficients[o3_coef_row, "Estimate"]
        se_value_o3   <- coefficients[o3_coef_row, "Std. Error"]
        p_value_o3    <- coefficients[o3_coef_row, "Pr(>|t|)"]

        ci_lower_o3 <- coef_value_o3 - 1.96 * se_value_o3
        ci_upper_o3 <- coef_value_o3 + 1.96 * se_value_o3

        coef_iqr_o3      <- coef_value_o3 * iqr_o3
        ci_lower_iqr_o3  <- ci_lower_o3 * iqr_o3
        ci_upper_iqr_o3  <- ci_upper_o3 * iqr_o3

        significance_o3 <- ""
        if (p_value_o3 < 0.001) {
          significance_o3 <- "***"
        } else if (p_value_o3 < 0.01) {
          significance_o3 <- "**"
        } else if (p_value_o3 < 0.05) {
          significance_o3 <- "*"
        }

        results$o3 <- list(
          exposure        = o3_name,
          adjustment      = adjustment_level,
          coefficient     = coef_value_o3,
          coefficient_iqr = coef_iqr_o3,
          ci_lower        = ci_lower_o3,
          ci_upper        = ci_upper_o3,
          ci_lower_iqr    = ci_lower_iqr_o3,
          ci_upper_iqr    = ci_upper_iqr_o3,
          p_value         = p_value_o3,
          significance    = significance_o3,
          buffer          = buffer_label,
          model_converged = TRUE
        )
      }

      return(results)

    } else {
      # Single-exposure model
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

      if (exposure_var == evi_name) {
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
        buffer          = buffer_label,
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
      buffer          = buffer_label,
      model_converged = FALSE,
      error_message   = e$message
    ))
  })
}

# =============================================================================
# 6. Run all models and collect results (1KM + 500m, with dual exposure)
# =============================================================================

adjustment_levels <- c("Unadjusted", "Age and sex adjusted", "Multivariable adjusted")

buffers_info <- list(
  "1KM" = list(
    evi_name = "EVI_1KM_mean4y",
    o3_name  = "O3_1KM_mean4y",
    iqr_evi  = iqr_evi_1km,
    iqr_o3   = iqr_o3_1km
  ),
  "500m" = list(
    evi_name = "EVI_500m_mean4y",
    o3_name  = "O3_500m_mean4y",
    iqr_evi  = iqr_evi_500,
    iqr_o3   = iqr_o3_500
  )
)

results <- list()
model_count <- 0
successful_models <- 0

for (buf in names(buffers_info)) {
  info <- buffers_info[[buf]]
  exposure_vars <- c(info$evi_name, info$o3_name)

  # Single-exposure models
  for (exp_var in exposure_vars) {
    for (adj_level in adjustment_levels) {
      model_count <- model_count + 1
      cat("Running single-exposure model", model_count, ": buffer", buf, ",",
          exp_var, " - ", adj_level, "\n")

      model_result <- run_lmm_model_safe(
        exposure_var   = exp_var,
        adjustment_level = adj_level,
        data           = data,
        evi_name       = info$evi_name,
        o3_name        = info$o3_name,
        iqr_evi        = info$iqr_evi,
        iqr_o3         = info$iqr_o3,
        buffer_label   = buf
      )

      if (!is.null(model_result) && model_result$model_converged) {
        res_name <- paste(buf, exp_var, adj_level, sep = "_")
        results[[res_name]] <- model_result
        successful_models <- successful_models + 1
      }
    }
  }

  # Dual-exposure models
  cat("Running dual-exposure model: buffer", buf, "EVI and O3 mutually adjusted\n")
  dual_exposure_result <- run_lmm_model_safe(
    exposure_var     = "dual",
    adjustment_level = "Dual exposure adjusted",
    data             = data,
    evi_name         = info$evi_name,
    o3_name          = info$o3_name,
    iqr_evi          = info$iqr_evi,
    iqr_o3           = info$iqr_o3,
    buffer_label     = buf
  )

  if (!is.null(dual_exposure_result)) {
    if (!is.null(dual_exposure_result$evi) &&
        dual_exposure_result$evi$model_converged) {
      results[[paste(buf, "EVI_Dual exposure adjusted", sep = "_")]] <- dual_exposure_result$evi
      successful_models <- successful_models + 1
    }
    if (!is.null(dual_exposure_result$o3) &&
        dual_exposure_result$o3$model_converged) {
      results[[paste(buf, "O3_Dual exposure adjusted", sep = "_")]] <- dual_exposure_result$o3
      successful_models <- successful_models + 1
    }
  }
}

total_expected <- (length(c("EVI","O3")) * length(adjustment_levels) + 2) *
  length(buffers_info)
cat("Successfully fitted coefficients:", successful_models, "/", total_expected, "\n")

# =============================================================================
# 7. Organize results and create tables
# =============================================================================

if (successful_models == 0) {
  stop("All model fitting failed; please check data quality and model specification")
}

results_list <- list()
for (result_name in names(results)) {
  result <- results[[result_name]]
  if (is.list(result) && !is.null(result$model_converged) &&
      result$model_converged) {
    results_list[[result_name]] <- as.data.frame(result)
  }
}

if (length(results_list) > 0) {
  results_df <- do.call(rbind, results_list)

  # Buffer column
  results_df$Buffer <- factor(results_df$buffer, levels = c("1KM","500m"))

  # Exposure factor (EVI / O3)
  results_df$Exposure_Factor <- ifelse(grepl("EVI", results_df$exposure), "EVI", "O3")

  # Model type (single / dual exposure)
  results_df$Model_Type <- ifelse(results_df$adjustment == "Dual exposure adjusted",
                            "Dual-exposure model", "Single-exposure model")
  results_df$Display_Factor <- ifelse(results_df$adjustment == "Dual exposure adjusted",
                            paste0(results_df$Exposure_Factor, " (dual)"),
                            results_df$Exposure_Factor)

  # Display string (per IQR effect + 95% CI)
  results_df$Display_Result <- ifelse(!is.na(results_df$coefficient_iqr),
                            sprintf("%.3f (%.3f-%.3f)%s",
                                    results_df$coefficient_iqr,
                                    results_df$ci_lower_iqr,
                                    results_df$ci_upper_iqr,
                                    results_df$significance),
                            "Model fitting failed")

  # Wide table
  wide_results <- results_df %>%
    select(Buffer, Exposure_Factor, adjustment, Display_Result) %>%
    tidyr::pivot_wider(names_from = adjustment, values_from = Display_Result)

} else {
  stop("No successfully fitted model results")
}

# =============================================================================
# 8. Create main results table (can output by Buffer separately)
# =============================================================================

if (nrow(wide_results) > 0) {
  main_table <- kable(wide_results,
                      col.names = c("Buffer", "Exposure Factor",
                                    "Unadjusted", "Age and sex adjusted",
                                    "Multivariable adjusted", "Dual exposure adjusted"),
                      align = c("c","l","c","c","c","c"),
                      caption = "Per-IQR LMM coefficient estimates (95% CI) of EVI and O3 exposure on BAA_PhenoAge across different buffers",
                      format = "html") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = FALSE) %>%
    add_header_above(c(" " = 2, "Adjustment Type" = 4)) %>%
    row_spec(0, bold = TRUE) %>%
    column_spec(2, bold = TRUE)

  print(main_table)
  save_kable(main_table, file = "output/tables/iqr_two_scales_results.html")
} else {
  warning("No successfully fitted model results to display")
}

# =============================================================================
# 9. Create visualizations (integrated 1KM + 500m)
# =============================================================================

if (nrow(results_df) > 0) {
  plot_data <- results_df %>%
    mutate(
      Adjustment_Type = factor(
        adjustment,
        levels = c("Unadjusted", "Age and sex adjusted", "Multivariable adjusted", "Dual exposure adjusted"),
        labels = c("Unadjusted", "Age and sex adjusted", "Multivariable adjusted", "Dual exposure adjusted")
      )
    )

  # ---------- Combined coefficient plot (points + error bars) ----------
  coefficient_plot <- ggplot(plot_data,
                             aes(x = Adjustment_Type,
                                 y = coefficient_iqr,
                                 color = Display_Factor,
                                 group = interaction(Display_Factor, Buffer))) +
    geom_point(position = position_dodge(width = 0.7), size = 2.6) +
    geom_errorbar(aes(ymin = ci_lower_iqr, ymax = ci_upper_iqr),
                  width = 0.2,
                  position = position_dodge(width = 0.7)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ Buffer, nrow = 1) +
    labs(title    = "Effect of EVI and O3 exposure on biological age acceleration: 1KM vs 500m comparison",
         subtitle = "Per-IQR coefficient estimates (95% CI), single- and dual-exposure models",
         x        = "Adjustment type",
         y        = "Coefficient estimate (per IQR increase)",
         color    = "Exposure factor") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.x     = element_text(angle = 45, hjust = 1)
    )

  print(coefficient_plot)
  ggsave("output/figures/iqr_two_scales_coefficient_plot.png",
         coefficient_plot, width = 12, height = 6, dpi = 300)

  # ---------- Combined group comparison plot (bar chart) ----------
  group_plot <- ggplot(plot_data,
                       aes(x = Exposure_Factor,
                           y = coefficient_iqr,
                           fill = Adjustment_Type)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = ci_lower_iqr, ymax = ci_upper_iqr),
                  position = position_dodge(width = 0.8),
                  width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ Buffer, nrow = 1) +
    labs(title    = "Exposure effect comparison: single vs dual-exposure models (1KM and 500m)",
         subtitle = "Per-IQR coefficient estimates (95% CI)",
         x        = "Exposure factor",
         y        = "Coefficient estimate (per IQR increase)",
         fill     = "Model type") +
    theme_minimal() +
    theme(legend.position = "bottom")

  print(group_plot)
  ggsave("output/figures/iqr_two_scales_group_comparison.png",
         group_plot, width = 12, height = 6, dpi = 300)
}

# =============================================================================
# 10. Save detailed results (merged 1KM + 500m)
# =============================================================================

write.csv(results_df,
          "output/tables/iqr_two_scales_detailed_results.csv",
          row.names = FALSE)

# =============================================================================
# 11. Output summary information
# =============================================================================

cat("\n")
cat("==================================================\n")
cat("Analysis complete! (1KM and 500m dual-scale, with dual-exposure models)\n")
cat("==================================================\n")
cat("Data information:\n")
cat("- Cleaned data dimensions:", nrow(data), "rows,", ncol(data), "columns\n")
cat("- 1KM EVI IQR:",  round(iqr_evi_1km, 3), "\n")
cat("- 1KM O3  IQR:",  round(iqr_o3_1km, 3), "\n")
cat("- 500m EVI IQR:", round(iqr_evi_500, 3), "\n")
cat("- 500m O3  IQR:", round(iqr_o3_500, 3), "\n")
cat("\n")
cat("Model fitting results:\n")
cat("- Single-exposure models (two buffers, total):",
    length(c("EVI","O3")) * length(adjustment_levels) * length(buffers_info),
    "\n")
cat("- Dual-exposure models:", length(buffers_info),
    "(each containing both EVI and O3 coefficients)\n")
cat("- Successfully fitted coefficients:", successful_models, "\n")
cat("\n")
cat("Results files saved:\n")
cat("- HTML table: ../output/tables/iqr_two_scales_results.html\n")
cat("- Coefficient plot: ../output/figures/iqr_two_scales_coefficient_plot.png\n")
cat("- Group comparison plot: ../output/figures/iqr_two_scales_group_comparison.png\n")
cat("- Detailed data: ../output/tables/iqr_two_scales_detailed_results.csv\n")
cat("\n")
cat("Significance notation:\n")
cat("* p < 0.05, ** p < 0.01, *** p < 0.001\n")
cat("==================================================\n")

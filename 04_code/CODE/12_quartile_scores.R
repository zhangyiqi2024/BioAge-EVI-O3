# =============================================================================
# Complete R code: Quartile scores and composite exposure score analysis
# =============================================================================

set.seed(20240101)

# Load required R packages
required_packages <- c("lme4", "lmerTest", "dplyr", "tidyr",
                       "kableExtra", "haven", "ggplot2", "scales")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# =============================================================================
# 1. Data loading and preprocessing
# =============================================================================

# Load data
data <- NULL
data_file_csv <- "DATA/analysis_ready/analysis_dataset.csv"
data_file_dta <- "DATA/analysis_ready/analysis_dataset.dta"

if (file.exists(data_file_csv)) {
  data <- read.csv(data_file_csv, stringsAsFactors = FALSE)
  cat("Successfully loaded CSV format data\n")
} else if (file.exists(data_file_dta)) {
  data <- as.data.frame(read_dta(data_file_dta))
  cat("Successfully loaded Stata format data\n")
}

if (is.null(data)) {
  stop("Failed to find data file. Please check: ", data_file_csv, " or ", data_file_dta)
}

cat("Original data dimensions:", dim(data), "\n")

# =============================================================================
# 2. Data cleaning
# =============================================================================

analysis_vars <- c("BAA_PhenoAge", "EVI_1KM_mean4y", "O3_1KM_mean4y",
                   "age", "sex", "id_num",
                   "education", "hh_consumption_total", "retirement", "pension_insurance",
                   "smoking_current", "alcohol_current",
                   "TMP", "PER", "SNOW", "DEM", "lat",
                   "bmi", "hypertension", "dyslipidemia", "year")

analysis_data <- data[intersect(analysis_vars, names(data))]

# Create log household consumption variable
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

# Process categorical variables
categorical_vars <- c("sex", "education", "retirement", "pension_insurance",
                      "smoking_current", "alcohol_current", "hypertension", "dyslipidemia", "year")

for (var in categorical_vars) {
  if (var %in% names(analysis_data)) {
    analysis_data[[var]] <- as.factor(analysis_data[[var]])
  }
}

# Remove rows with any invalid values
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

data <- analysis_data_clean

# =============================================================================
# 3. Climate variable principal component analysis
# =============================================================================

cat("=== Climate Variable PCA ===\n")

climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")
climate_data <- data[climate_vars]

if (any(is.na(climate_data))) {
  cat("Climate variables have missing values, will be removed\n")
  complete_cases <- complete.cases(climate_data)
  climate_data <- climate_data[complete_cases, ]
  data <- data[complete_cases, ]
}

pca_result <- prcomp(climate_data, scale. = TRUE, center = TRUE)
print(summary(pca_result))

climate_pcs <- as.data.frame(pca_result$x[, 1:5])
colnames(climate_pcs) <- paste0("climate_pc", 1:5)
data <- cbind(data, climate_pcs)

cat("Climate PCA complete, 5 PCs added to data\n")

# =============================================================================
# 4. Calculate quartile scores for exposure variables
# =============================================================================

calculate_quartile_scores <- function(x, reverse = FALSE) {
  # Calculate quartiles
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

  # Determine score direction based on reverse parameter
  if (reverse) {
    # Reverse scoring: lower values get higher scores
    scores <- ifelse(x <= q[1], 3,
                     ifelse(x <= q[2], 2,
                            ifelse(x <= q[3], 1, 0)))
  } else {
    # Forward scoring: higher values get higher scores
    scores <- ifelse(x <= q[1], 0,
                     ifelse(x <= q[2], 1,
                            ifelse(x <= q[3], 2, 3)))
  }
  return(scores)
}

# Calculate quartile scores
# EVI: positive indicator, higher values get higher scores
data$EVI_score <- calculate_quartile_scores(data$EVI_1KM_mean4y, reverse = FALSE)
# O3: negative indicator, lower values get higher scores
data$O3_score <- calculate_quartile_scores(data$O3_1KM_mean4y, reverse = TRUE)

# Calculate mixed exposure score (0-6 points)
data$mixed_exposure_score <- data$EVI_score + data$O3_score

# Output score distributions
cat("EVI score distribution:\n")
print(table(data$EVI_score))
cat("O3 score distribution:\n")
print(table(data$O3_score))
cat("Mixed exposure score distribution:\n")
print(table(data$mixed_exposure_score))

# =============================================================================
# 5. Linear mixed model function (using score variables)
# =============================================================================

run_lmm_model_scores <- function(exposure_type, adjustment_level, data) {
  tryCatch({
    if (exposure_type == "EVI_score") {
      exposure_var <- "EVI_score"
    } else if (exposure_type == "O3_score") {
      exposure_var <- "O3_score"
    } else if (exposure_type == "mixed_score") {
      exposure_var <- "mixed_exposure_score"
    } else if (exposure_type == "dual_scores") {
      # Dual-exposure model: includes both EVI and O3 scores
      exposure_var <- "EVI_score + O3_score"
    }

    if (adjustment_level == "Unadjusted") {
      formula_str <- paste("BAA_PhenoAge ~", exposure_var, "+ (1 | id_num)")
    } else if (adjustment_level == "Age and sex adjusted") {
      formula_str <- paste("BAA_PhenoAge ~", exposure_var, "+ age + sex + (1 | id_num)")
    } else if (adjustment_level == "Multivariable adjusted") {
      formula_str <- paste("BAA_PhenoAge ~", exposure_var,
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

    if (exposure_type == "dual_scores") {
      results <- list()

      # EVI score
      evi_coef_row <- which(rownames(coefficients) == "EVI_score")
      if (length(evi_coef_row) > 0) {
        coef_value_evi <- coefficients[evi_coef_row, "Estimate"]
        se_value_evi   <- coefficients[evi_coef_row, "Std. Error"]
        p_value_evi    <- coefficients[evi_coef_row, "Pr(>|t|)"]

        ci_lower_evi <- coef_value_evi - 1.96 * se_value_evi
        ci_upper_evi <- coef_value_evi + 1.96 * se_value_evi

        significance_evi <- ""
        if (p_value_evi < 0.001) {
          significance_evi <- "***"
        } else if (p_value_evi < 0.01) {
          significance_evi <- "**"
        } else if (p_value_evi < 0.05) {
          significance_evi <- "*"
        }

        results$evi <- list(
          exposure        = "EVI_score",
          adjustment      = adjustment_level,
          coefficient     = coef_value_evi,
          ci_lower        = ci_lower_evi,
          ci_upper        = ci_upper_evi,
          p_value         = p_value_evi,
          significance    = significance_evi,
          model_converged = TRUE
        )
      }

      # O3 score
      o3_coef_row <- which(rownames(coefficients) == "O3_score")
      if (length(o3_coef_row) > 0) {
        coef_value_o3 <- coefficients[o3_coef_row, "Estimate"]
        se_value_o3   <- coefficients[o3_coef_row, "Std. Error"]
        p_value_o3    <- coefficients[o3_coef_row, "Pr(>|t|)"]

        ci_lower_o3 <- coef_value_o3 - 1.96 * se_value_o3
        ci_upper_o3 <- coef_value_o3 + 1.96 * se_value_o3

        significance_o3 <- ""
        if (p_value_o3 < 0.001) {
          significance_o3 <- "***"
        } else if (p_value_o3 < 0.01) {
          significance_o3 <- "**"
        } else if (p_value_o3 < 0.05) {
          significance_o3 <- "*"
        }

        results$o3 <- list(
          exposure        = "O3_score",
          adjustment      = adjustment_level,
          coefficient     = coef_value_o3,
          ci_lower        = ci_lower_o3,
          ci_upper        = ci_upper_o3,
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
        ci_lower        = ci_lower,
        ci_upper        = ci_upper,
        p_value         = p_value,
        significance    = significance,
        model_converged = TRUE
      ))
    }
  }, error = function(e) {
    warning("Model fitting failed: ", e$message)
    return(list(
      exposure        = exposure_type,
      adjustment      = adjustment_level,
      coefficient     = NA,
      ci_lower        = NA,
      ci_upper        = NA,
      p_value         = NA,
      significance    = "",
      model_converged = FALSE,
      error_message   = e$message
    ))
  })
}

# =============================================================================
# 6. Run all models and collect results (fixed version)
# =============================================================================

exposure_types <- c("EVI_score", "O3_score", "mixed_score", "dual_scores")
adjustment_levels <- c("Unadjusted", "Age and sex adjusted", "Multivariable adjusted")

results <- list()
model_count <- 0
successful_models <- 0

for (exp_type in exposure_types) {
  for (adj_level in adjustment_levels) {
    model_count <- model_count + 1
    cat("Running model", model_count, ": ", exp_type, " - ", adj_level, "\n")

    model_result <- run_lmm_model_scores(exp_type, adj_level, data)

    # Fix: more robust error checking
    if (!is.null(model_result)) {
      if (exp_type == "dual_scores") {
        # Special handling for dual-exposure models
        if (!is.null(model_result$evi) &&
            is.list(model_result$evi) &&
            !is.null(model_result$evi$model_converged) &&
            model_result$evi$model_converged) {
          results[[paste("EVI_score", adj_level, "dual", sep = "_")]] <- model_result$evi
          successful_models <- successful_models + 1
          cat("  EVI score model succeeded\n")
        } else {
          cat("  EVI score model failed\n")
        }

        if (!is.null(model_result$o3) &&
            is.list(model_result$o3) &&
            !is.null(model_result$o3$model_converged) &&
            model_result$o3$model_converged) {
          results[[paste("O3_score", adj_level, "dual", sep = "_")]] <- model_result$o3
          successful_models <- successful_models + 1
          cat("  O3 score model succeeded\n")
        } else {
          cat("  O3 score model failed\n")
        }
      } else {
        # Single-exposure model handling
        if (is.list(model_result) &&
            !is.null(model_result$model_converged) &&
            model_result$model_converged) {
          results[[paste(exp_type, adj_level, sep = "_")]] <- model_result
          successful_models <- successful_models + 1
          cat("  Model succeeded\n")
        } else {
          cat("  Model failed\n")
        }
      }
    } else {
      cat("  Model returned NULL\n")
    }
  }
}

cat("Number of successfully fitted models:", successful_models, "/",
    (length(exposure_types) * length(adjustment_levels) + 2), "\n")

# =============================================================================
# 7. Organize results and create tables
# =============================================================================

if (successful_models == 0) {
  stop("All models failed to fit. Please check data quality and model specifications.")
}

results_list <- list()
for (result_name in names(results)) {
  result <- results[[result_name]]
  if (is.list(result) && !is.null(result$model_converged) && result$model_converged) {
    results_list[[result_name]] <- as.data.frame(result)
  }
}

if (length(results_list) == 0) {
  stop("No successfully fitted model results available")
}

results_df <- do.call(rbind, results_list)

# Set exposure factor labels
results_df$exposure_factor <- ifelse(grepl("EVI_score", results_df$exposure), "EVI Score",
                              ifelse(grepl("O3_score", results_df$exposure), "O₃ Score",
                                     ifelse(grepl("mixed", results_df$exposure), "Mixed Exposure Score", "Unknown")))

results_df$model_type <- ifelse(grepl("dual", rownames(results_df)), "Dual-exposure", "Single-exposure")

# Calculate SE (back-calculated from CI)
results_df <- results_df %>%
  mutate(
    SE = (ci_upper - ci_lower) / (2 * 1.96),
    Exposure_label = case_when(
      exposure_factor == "EVI Score" ~ "EVI_score",
      exposure_factor == "O₃ Score" ~ "O3_score",
      exposure_factor == "Mixed Exposure Score" ~ "mixed_score",
      TRUE ~ exposure
    ),
    ModelType_en = ifelse(model_type == "Dual-exposure", "Dual-scores", "Single-score"),
    Buffer = "1KM",
    Adjustment_en = adjustment
  )

# Generate display column for main HTML table
results_df <- results_df %>%
  mutate(
    Display = ifelse(
      !is.na(coefficient),
      sprintf("%.3f (%.3f, %.3f)%s",
              coefficient,
              ci_lower,
              ci_upper,
              significance),
      "Model failed"
    )
  )

# =============================================================================
# 8. Create main results HTML table
# =============================================================================

# Reorganize results table
wide_results <- results_df %>%
  select(Exposure = exposure_factor,
         Adjustment = Adjustment_en,
         Display) %>%
  tidyr::pivot_wider(names_from = Adjustment, values_from = Display)

desired_columns <- c("Exposure", "Unadjusted", "Age and sex adjusted", "Multivariable adjusted")
existing_columns <- intersect(desired_columns, names(wide_results))
wide_results <- wide_results[, existing_columns, drop = FALSE]

if (nrow(wide_results) > 0) {
  main_table <- kable(
    wide_results,
    col.names = c("Exposure", "Unadjusted", "Age and sex adjusted", "Multivariable adjusted"),
    align = c("l", "c", "c", "c"),
    caption = "Associations of EVI score, O₃ score, and mixed exposure score with phenotypic age acceleration (per 1-point increase in score)",
    format = "html"
  ) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = FALSE) %>%
    add_header_above(c(" " = 1, "Adjustment set" = 3)) %>%
    row_spec(0, bold = TRUE) %>%
    column_spec(1, bold = TRUE)

  print(main_table)
  save_kable(main_table, file = "output/tables/quartile_score_analysis.html")
} else {
  warning("No successfully fitted model results to display")
}

# =============================================================================
# 9. Create visualizations (high-impact journal standard - English labels)
# =============================================================================

if (nrow(results_df) > 0) {

  # Create professional English labels
  plot_data <- results_df %>%
    mutate(
      Adjustment = factor(
        Adjustment_en,
        levels = c("Unadjusted", "Age and sex adjusted", "Multivariable adjusted")
      ),
      # Translate exposure factor labels to English
      Exposure_en = case_when(
        exposure_factor == "EVI Score" ~ "EVI Score",
        exposure_factor == "O₃ Score" ~ "O₃ Score",
        exposure_factor == "Mixed Exposure Score" ~ "Mixed Exposure Score",
        TRUE ~ exposure_factor
      ),
      Model_Type = model_type,
      # Professional scientific labels (English)
      Exposure_Label = case_when(
        model_type == "Dual-exposure" & exposure_factor == "EVI Score" ~ "EVI Score (Dual)",
        model_type == "Dual-exposure" & exposure_factor == "O₃ Score" ~ "O₃ Score (Dual)",
        model_type == "Single-exposure" & exposure_factor == "EVI Score" ~ "EVI Score (Single)",
        model_type == "Single-exposure" & exposure_factor == "O₃ Score" ~ "O₃ Score (Single)",
        model_type == "Single-exposure" & exposure_factor == "Mixed Exposure Score" ~ "Mixed Exposure Score",
        TRUE ~ exposure_factor
      ),
      # Model type English labels
      Model_Type_en = case_when(
        model_type == "Single-exposure" ~ "Single-Exposure",
        model_type == "Dual-exposure" ~ "Dual-Exposure",
        TRUE ~ model_type
      )
    )

  # Figure 1: Forest plot - English labels
  coefficient_plot <- ggplot(
    plot_data,
    aes(x = Exposure_Label, y = coefficient,
        color = Model_Type_en, shape = Model_Type_en)
  ) +
    geom_point(position = position_dodge(width = 0.7), size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                  width = 0.2,
                  position = position_dodge(width = 0.7)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.7) +
    facet_wrap(~ Adjustment, ncol = 1) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
    scale_color_manual(values = c("Single-Exposure" = "#1f77b4",
                                  "Dual-Exposure" = "#ff7f0e")) +
    scale_shape_manual(values = c("Single-Exposure" = 16,
                                  "Dual-Exposure" = 17)) +
    labs(
      title = "Associations of Environmental Exposure Scores with Phenotypic Age Acceleration",
      subtitle = "Linear Mixed-Effects Models: Single- vs. Dual-Exposure Approaches",
      x = "Exposure Score Type",
      y = "β Coefficient (per 1-Score Increase)",
      color = "Model Type",
      shape = "Model Type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12)
    ) +
    coord_flip()

  print(coefficient_plot)
  ggsave("output/figures/quartile_scores_forest_plot.png",
         coefficient_plot, width = 10, height = 8, dpi = 300)

  # Figure 2: Grouped bar chart - English labels
  group_plot <- ggplot(
    plot_data,
    aes(x = Adjustment, y = coefficient, fill = Model_Type_en)
  ) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                  position = position_dodge(width = 0.8),
                  width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ Exposure_en, scales = "free_y") +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
    scale_fill_manual(values = c("Single-Exposure" = "#1f77b4",
                                 "Dual-Exposure" = "#ff7f0e")) +
    labs(
      title = "Comparative Effects of Environmental Exposures on Phenotypic Aging",
      subtitle = "Quartile-Based Exposure Scores Across Adjustment Models",
      x = "Adjustment Set",
      y = "β Coefficient (per 1-Score Increase)",
      fill = "Model Type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey95", color = NA),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12)
    )

  print(group_plot)
  ggsave("output/figures/quartile_scores_grouped_bar.png",
         group_plot, width = 12, height = 6, dpi = 300)

  # Figure 3: Model comparison plot - English labels
  comparison_data <- plot_data %>%
    filter(Exposure_en %in% c("EVI Score", "O₃ Score"))

  if (nrow(comparison_data) > 0) {
    comparison_plot <- ggplot(
      comparison_data,
      aes(x = Adjustment, y = coefficient, color = Model_Type_en, group = Model_Type_en)
    ) +
      geom_line(position = position_dodge(width = 0.3), alpha = 0.6) +
      geom_point(position = position_dodge(width = 0.3), size = 3) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                    width = 0.2,
                    position = position_dodge(width = 0.3)) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      facet_wrap(~ Exposure_en, ncol = 2) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
      scale_color_manual(values = c("Single-Exposure" = "#1f77b4",
                                    "Dual-Exposure" = "#ff7f0e")) +
      labs(
        title = "Modeling Approach Comparison: Single vs. Dual Environmental Exposures",
        subtitle = "Coefficient Trajectories Across Adjustment Sets for EVI and O₃ Scores",
        x = "Adjustment Set",
        y = "β Coefficient (per 1-Score Increase)",
        color = "Model Type"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey95", color = NA),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 12)
      )

    print(comparison_plot)
    ggsave("output/figures/quartile_scores_model_comparison.png",
           comparison_plot, width = 10, height = 6, dpi = 300)
  }

  # Figure 4: Mixed exposure score dedicated analysis (English labels)
  mixed_plot_data <- plot_data %>%
    filter(Exposure_en == "Mixed Exposure Score")

  if (nrow(mixed_plot_data) > 0) {
    mixed_plot <- ggplot(
      mixed_plot_data,
      aes(x = Adjustment, y = coefficient)
    ) +
      geom_col(fill = "#2ca02c", alpha = 0.8, width = 0.6) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                    width = 0.2, color = "darkgreen") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_text(aes(label = sprintf("%.3f", coefficient)),
                vjust = -0.5, size = 3.5, fontface = "bold") +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
      labs(
        title = "Combined Environmental Exposure Score and Phenotypic Aging",
        subtitle = "Mixed Exposure Score (EVI + O₃ Quartiles) Association with Age Acceleration",
        x = "Adjustment Set",
        y = "β Coefficient (per 1-Score Increase)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 12),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()
      )

    print(mixed_plot)
    ggsave("output/figures/quartile_scores_mixed_exposure.png",
           mixed_plot, width = 8, height = 6, dpi = 300)
  }
}

# =============================================================================
# 10. Save detailed results
# =============================================================================

results_for_csv <- results_df %>%
  mutate(
    Exposure = Exposure_label,
    `Model Type` = ModelType_en,
    Buffer = Buffer,
    Adjustment_csv = dplyr::case_when(
      adjustment == "Unadjusted" ~ "Unadjusted",
      adjustment == "Age and sex adjusted" ~ "Age and sex",
      adjustment == "Multivariable adjusted" ~ "Multivariable",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(
    factor(Exposure, levels = c("EVI_score", "O3_score", "mixed_score")),
    factor(`Model Type`, levels = c("Single-score", "Dual-scores")),
    factor(Adjustment_csv, levels = c("Unadjusted", "Age and sex", "Multivariable"))
  ) %>%
  transmute(
    Exposure,
    `Model Type`,
    Buffer,
    Adjustment = Adjustment_csv,
    `Coefficient (SE)` = sprintf("%.3f (%.3f)", coefficient, SE),
    p_value = sprintf("%.3f", p_value),
    `95% CI` = sprintf("[%.3f, %.3f]", ci_lower, ci_upper)
  )

write.csv(results_for_csv, "output/tables/quartile_scores_detailed_results.csv", row.names = FALSE)

# =============================================================================
# 11. Output summary information
# =============================================================================

cat("\n")
cat("==================================================\n")
cat("Quartile score analysis complete!\n")
cat("==================================================\n")
cat("Score construction information:\n")
cat("- EVI score: 0-3 points (higher value = higher score)\n")
cat("- O₃ score: 0-3 points (lower value = higher score)\n")
cat("- Mixed exposure score: 0-6 points (EVI score + O₃ score)\n")
cat("\n")
cat("Score distributions:\n")
cat("- EVI score: ", paste(names(table(data$EVI_score)), collapse = ", "), "\n")
cat("- O₃ score: ", paste(names(table(data$O3_score)), collapse = ", "), "\n")
cat("- Mixed exposure score: ", paste(names(table(data$mixed_exposure_score)), collapse = ", "), "\n")
cat("\n")
cat("Data information:\n")
cat("- Cleaned data dimensions:", nrow(data), "rows,", ncol(data), "columns\n")
cat("\n")
cat("Model fitting results:\n")
cat("- Number of successfully fitted models:", successful_models, "\n")
cat("\n")
cat("Results files saved:\n")
cat("- HTML table: ../output/tables/quartile_score_analysis.html\n")
cat("- Forest plot: ../output/figures/quartile_scores_forest_plot.png\n")
cat("- Grouped bar chart: ../output/figures/quartile_scores_grouped_bar.png\n")
cat("- Model comparison: ../output/figures/quartile_scores_model_comparison.png\n")
cat("- Mixed exposure: ../output/figures/quartile_scores_mixed_exposure.png\n")
cat("- Detailed CSV: ../output/tables/quartile_scores_detailed_results.csv\n")
cat("\n")
cat("Significance notation:\n")
cat("* p < 0.05, ** p < 0.01, *** p < 0.001\n")
cat("==================================================\n")

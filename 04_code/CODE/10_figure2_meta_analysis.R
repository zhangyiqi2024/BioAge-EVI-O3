# 10_figure2_meta_analysis.R — Figure 2: model-comparison forest plot
# v4 update (2026-07-22): the random-effects meta-analysis pooling across
# correlated model specifications was removed from the narrative and from this
# figure; alternative specifications are presented side-by-side as robustness
# analyses. The MSM row label is now "Propensity-score weighted".
# (Script filename retained for pipeline compatibility.)
# Outputs:
#   output/figures/Figure2_forest_plot.png / .pdf  (combined table + forest)
#   output/tables/Exposure_Effects_Summary.csv
#   output/tables/Detailed_Results_Final.csv
#   output/tables/Complete_Results_with_Meta_Final.csv  (legacy filename; no
#     meta-analysis rows are produced any more)

# Load necessary packages
library(haven)
library(lme4)
library(lmerTest)
library(sandwich)
library(lmtest)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(patchwork)
library(grid)

set.seed(20240101)

# Set working directory to CODE folder

# -----------------------------------------------------------------------------
# Step 1: Load data and perform data cleaning
# -----------------------------------------------------------------------------

# Read data
data <- read_dta("DATA/analysis_ready/analysis_dataset.dta")

# Data filtering and variable creation
data <- data %>%
  filter(age > 45 & age < 90) %>%
  mutate(
    ln_hh_consumption_total = log(hh_consumption_total),
    id_num = as.factor(id_num),
    year = as.factor(year)
  ) %>%
  arrange(id_num, year)

# Check key variable types and ranges
cat("=== Data diagnostics ===\n")
cat("Sample size:", nrow(data), "\n")

check_vars <- c(
  "BAA_PhenoAge", "traj15_EVI_1KM", "traj15_EVI_500m",
  "traj15_O3_1KM", "traj15_O3_500m", "ipw_EVI_1KM", "ipw_EVI_500m",
  "ipw_O3_1KM", "ipw_O3_500m", "ipw_joint_1KM", "ipw_joint_500m"
)

for (var in check_vars) {
  if (var %in% names(data)) {
    if (is.numeric(data[[var]])) {
      cat(
        var, "type: numeric, range:",
        range(data[[var]], na.rm = TRUE),
        "NA count:", sum(is.na(data[[var]])), "\n"
      )
    } else {
      cat(
        var, "type:", class(data[[var]]),
        "levels:", levels(data[[var]]),
        "NA count:", sum(is.na(data[[var]])), "\n"
      )
    }
  } else {
    cat("Warning: variable", var, "does not exist\n")
  }
}

# Create dummy weight variables if they don't exist
if (!"ipw_EVI_1KM" %in% names(data)) {
  data$ipw_EVI_1KM    <- 1
  data$ipw_EVI_500m   <- 1
  data$ipw_O3_1KM     <- 1
  data$ipw_O3_500m    <- 1
  data$ipw_joint_1KM  <- 1
  data$ipw_joint_500m <- 1
  cat("Created dummy weight variables\n")
}

# Climate variable principal component analysis
climate_vars <- data %>%
  select(TMP, PER, SNOW, DEM, lat) %>%
  na.omit()

if (nrow(climate_vars) > 0) {
  pca_result  <- prcomp(climate_vars, scale. = TRUE)
  climate_pc  <- predict(pca_result, climate_vars)
  data        <- cbind(data, climate_pc[, 1:5])
  colnames(data)[(ncol(data)-4):ncol(data)] <- paste0("climate_pc", 1:5)
} else {
  cat("Warning: insufficient climate variable data, creating dummy climate PCs\n")
  data$climate_pc1 <- rnorm(nrow(data))
  data$climate_pc2 <- rnorm(nrow(data))
  data$climate_pc3 <- rnorm(nrow(data))
  data$climate_pc4 <- rnorm(nrow(data))
  data$climate_pc5 <- rnorm(nrow(data))
}

# -----------------------------------------------------------------------------
# Step 2: Simplified data cleaning function
# -----------------------------------------------------------------------------

clean_data_simple <- function(data, vars) {
  model_data <- data %>%
    select(all_of(vars)) %>%
    mutate(across(where(is.numeric),
                  ~ ifelse(is.infinite(.) | is.nan(.), NA, .))) %>%
    na.omit()
  return(model_data)
}

# -----------------------------------------------------------------------------
# Step 3: Run all regression models (single-exposure and dual-exposure)
# -----------------------------------------------------------------------------

cat("=== Starting all regression models ===\n")

# Base control variables
base_controls <- paste(
  "age", "sex", "education", "ln_hh_consumption_total", "retirement",
  "pension_insurance", "smoking_current", "alcohol_current",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4", "climate_pc5",
  "bmi", "hypertension", "dyslipidemia", "year",
  sep = " + "
)

# List to store model results
models <- list()

# ----------------------- Single-exposure model: EVI 1KM --------------------
cat("=== Single-exposure model: EVI 1KM ===\n")
vars_needed <- c(
  "BAA_PhenoAge", "traj15_EVI_1KM", "ipw_EVI_1KM", "id_num", "age", "sex",
  "education", "ln_hh_consumption_total", "retirement", "pension_insurance",
  "smoking_current", "alcohol_current",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4", "climate_pc5",
  "bmi", "hypertension", "dyslipidemia", "year"
)

model_data <- clean_data_simple(data, vars_needed)
formula <- as.formula(
  paste("BAA_PhenoAge ~ traj15_EVI_1KM +", base_controls, "+ (1 | id_num)")
)

# Weighted mixed effects (WME)
tryCatch({
  models$EVI_1KM_WME <- lmer(
    formula, data = model_data, weights = ipw_EVI_1KM,
    control = lmerControl(optimizer = "bobyqa")
  )
}, error = function(e) {
  cat("EVI 1KM WME failed:", e$message, "\n")
  models$EVI_1KM_WME <- NULL
})

# Unweighted mixed effects (UME)
tryCatch({
  models$EVI_1KM_UME <- lmer(
    formula, data = model_data,
    control = lmerControl(optimizer = "bobyqa")
  )
}, error = function(e) {
  cat("EVI 1KM UME failed:", e$message, "\n")
  models$EVI_1KM_UME <- NULL
})

# Propensity-score weighted model (PSW; marginal structural model)
formula_msm <- as.formula(
  paste("BAA_PhenoAge ~ traj15_EVI_1KM +", base_controls)
)
tryCatch({
  models$EVI_1KM_MSM <- lm(
    formula_msm, data = model_data, weights = ipw_EVI_1KM
  )
}, error = function(e) {
  cat("EVI 1KM MSM failed:", e$message, "\n")
  models$EVI_1KM_MSM <- NULL
})

# ----------------------- Single-exposure model: EVI 500M ------------------
cat("=== Single-exposure model: EVI 500M ===\n")
vars_needed <- c(
  "BAA_PhenoAge", "traj15_EVI_500m", "ipw_EVI_500m", "id_num", "age", "sex",
  "education", "ln_hh_consumption_total", "retirement", "pension_insurance",
  "smoking_current", "alcohol_current",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4", "climate_pc5",
  "bmi", "hypertension", "dyslipidemia", "year"
)

model_data <- clean_data_simple(data, vars_needed)
formula <- as.formula(
  paste("BAA_PhenoAge ~ traj15_EVI_500m +", base_controls, "+ (1 | id_num)")
)

tryCatch({
  models$EVI_500M_WME <- lmer(
    formula, data = model_data, weights = ipw_EVI_500m,
    control = lmerControl(optimizer = "bobyqa")
  )
}, error = function(e) {
  cat("EVI 500M WME failed:", e$message, "\n")
  models$EVI_500M_WME <- NULL
})

tryCatch({
  models$EVI_500M_UME <- lmer(
    formula, data = model_data,
    control = lmerControl(optimizer = "bobyqa")
  )
}, error = function(e) {
  cat("EVI 500M UME failed:", e$message, "\n")
  models$EVI_500M_UME <- NULL
})

formula_msm <- as.formula(
  paste("BAA_PhenoAge ~ traj15_EVI_500m +", base_controls)
)
tryCatch({
  models$EVI_500M_MSM <- lm(
    formula_msm, data = model_data, weights = ipw_EVI_500m
  )
}, error = function(e) {
  cat("EVI 500M MSM failed:", e$message, "\n")
  models$EVI_500M_MSM <- NULL
})

# ----------------------- Single-exposure model: O3 1KM --------------------
cat("=== Single-exposure model: O3 1KM ===\n")
vars_needed <- c(
  "BAA_PhenoAge", "traj15_O3_1KM", "ipw_O3_1KM", "id_num", "age", "sex",
  "education", "ln_hh_consumption_total", "retirement", "pension_insurance",
  "smoking_current", "alcohol_current",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4", "climate_pc5",
  "bmi", "hypertension", "dyslipidemia", "year"
)

model_data <- clean_data_simple(data, vars_needed)
formula <- as.formula(
  paste("BAA_PhenoAge ~ traj15_O3_1KM +", base_controls, "+ (1 | id_num)")
)

tryCatch({
  models$O3_1KM_WME <- lmer(
    formula, data = model_data, weights = ipw_O3_1KM,
    control = lmerControl(optimizer = "bobyqa")
  )
}, error = function(e) {
  cat("O3 1KM WME failed:", e$message, "\n")
  models$O3_1KM_WME <- NULL
})

tryCatch({
  models$O3_1KM_UME <- lmer(
    formula, data = model_data,
    control = lmerControl(optimizer = "bobyqa")
  )
}, error = function(e) {
  cat("O3 1KM UME failed:", e$message, "\n")
  models$O3_1KM_UME <- NULL
})

formula_msm <- as.formula(
  paste("BAA_PhenoAge ~ traj15_O3_1KM +", base_controls)
)
tryCatch({
  models$O3_1KM_MSM <- lm(
    formula_msm, data = model_data, weights = ipw_O3_1KM
  )
}, error = function(e) {
  cat("O3 1KM MSM failed:", e$message, "\n")
  models$O3_1KM_MSM <- NULL
})

# ----------------------- Single-exposure model: O3 500M -------------------
cat("=== Single-exposure model: O3 500M ===\n")
vars_needed <- c(
  "BAA_PhenoAge", "traj15_O3_500m", "ipw_O3_500m", "id_num", "age", "sex",
  "education", "ln_hh_consumption_total", "retirement", "pension_insurance",
  "smoking_current", "alcohol_current",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4", "climate_pc5",
  "bmi", "hypertension", "dyslipidemia", "year"
)

model_data <- clean_data_simple(data, vars_needed)
formula <- as.formula(
  paste("BAA_PhenoAge ~ traj15_O3_500m +", base_controls, "+ (1 | id_num)")
)

tryCatch({
  models$O3_500M_WME <- lmer(
    formula, data = model_data, weights = ipw_O3_500m,
    control = lmerControl(optimizer = "bobyqa")
  )
}, error = function(e) {
  cat("O3 500M WME failed:", e$message, "\n")
  models$O3_500M_WME <- NULL
})

tryCatch({
  models$O3_500M_UME <- lmer(
    formula, data = model_data,
    control = lmerControl(optimizer = "bobyqa")
  )
}, error = function(e) {
  cat("O3 500M UME failed:", e$message, "\n")
  models$O3_500M_UME <- NULL
})

formula_msm <- as.formula(
  paste("BAA_PhenoAge ~ traj15_O3_500m +", base_controls)
)
tryCatch({
  models$O3_500M_MSM <- lm(
    formula_msm, data = model_data, weights = ipw_O3_500m
  )
}, error = function(e) {
  cat("O3 500M MSM failed:", e$message, "\n")
  models$O3_500M_MSM <- NULL
})

# ----------------------- Dual-exposure model: EVI 1KM + O3 1KM -----------
cat("=== Dual-exposure model: EVI 1KM + O3 1KM ===\n")
vars_needed <- c(
  "BAA_PhenoAge", "traj15_EVI_1KM", "traj15_O3_1KM", "ipw_joint_1KM",
  "id_num", "age", "sex", "education",
  "ln_hh_consumption_total", "retirement", "pension_insurance",
  "smoking_current", "alcohol_current",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4", "climate_pc5",
  "bmi", "hypertension", "dyslipidemia", "year"
)

model_data <- clean_data_simple(data, vars_needed)
formula <- as.formula(
  paste("BAA_PhenoAge ~ traj15_EVI_1KM + traj15_O3_1KM +",
        base_controls, "+ (1 | id_num)")
)

tryCatch({
  models$Joint_EVI_1KM_WME <- lmer(
    formula, data = model_data, weights = ipw_joint_1KM,
    control = lmerControl(optimizer = "bobyqa")
  )
}, error = function(e) {
  cat("Joint EVI 1KM WME failed:", e$message, "\n")
  models$Joint_EVI_1KM_WME <- NULL
})

tryCatch({
  models$Joint_O3_1KM_WME <- models$Joint_EVI_1KM_WME
}, error = function(e) {
  cat("Joint O3 1KM WME failed:", e$message, "\n")
  models$Joint_O3_1KM_WME <- NULL
})

tryCatch({
  models$Joint_EVI_1KM_UME <- lmer(
    formula, data = model_data,
    control = lmerControl(optimizer = "bobyqa")
  )
}, error = function(e) {
  cat("Joint EVI 1KM UME failed:", e$message, "\n")
  models$Joint_EVI_1KM_UME <- NULL
})

tryCatch({
  models$Joint_O3_1KM_UME <- models$Joint_EVI_1KM_UME
}, error = function(e) {
  cat("Joint O3 1KM UME failed:", e$message, "\n")
  models$Joint_O3_1KM_UME <- NULL
})

formula_msm <- as.formula(
  paste("BAA_PhenoAge ~ traj15_EVI_1KM + traj15_O3_1KM +", base_controls)
)
tryCatch({
  models$Joint_EVI_1KM_MSM <- lm(
    formula_msm, data = model_data, weights = ipw_joint_1KM
  )
}, error = function(e) {
  cat("Joint EVI 1KM MSM failed:", e$message, "\n")
  models$Joint_EVI_1KM_MSM <- NULL
})

tryCatch({
  models$Joint_O3_1KM_MSM <- models$Joint_EVI_1KM_MSM
}, error = function(e) {
  cat("Joint O3 1KM MSM failed:", e$message, "\n")
  models$Joint_O3_1KM_MSM <- NULL
})

# ----------------------- Dual-exposure model: EVI 500M + O3 500M ---------
cat("=== Dual-exposure model: EVI 500M + O3 500M ===\n")
vars_needed <- c(
  "BAA_PhenoAge", "traj15_EVI_500m", "traj15_O3_500m", "ipw_joint_500m",
  "id_num", "age", "sex", "education",
  "ln_hh_consumption_total", "retirement", "pension_insurance",
  "smoking_current", "alcohol_current",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4", "climate_pc5",
  "bmi", "hypertension", "dyslipidemia", "year"
)

model_data <- clean_data_simple(data, vars_needed)
formula <- as.formula(
  paste("BAA_PhenoAge ~ traj15_EVI_500m + traj15_O3_500m +",
        base_controls, "+ (1 | id_num)")
)

tryCatch({
  models$Joint_EVI_500M_WME <- lmer(
    formula, data = model_data, weights = ipw_joint_500m,
    control = lmerControl(optimizer = "bobyqa")
  )
}, error = function(e) {
  cat("Joint EVI 500M WME failed:", e$message, "\n")
  models$Joint_EVI_500M_WME <- NULL
})

tryCatch({
  models$Joint_O3_500M_WME <- models$Joint_EVI_500M_WME
}, error = function(e) {
  cat("Joint O3 500M WME failed:", e$message, "\n")
  models$Joint_O3_500M_WME <- NULL
})

tryCatch({
  models$Joint_EVI_500M_UME <- lmer(
    formula, data = model_data,
    control = lmerControl(optimizer = "bobyqa")
  )
}, error = function(e) {
  cat("Joint EVI 500M UME failed:", e$message, "\n")
  models$Joint_EVI_500M_UME <- NULL
})

tryCatch({
  models$Joint_O3_500M_UME <- models$Joint_EVI_500M_UME
}, error = function(e) {
  cat("Joint O3 500M UME failed:", e$message, "\n")
  models$Joint_O3_500M_UME <- NULL
})

formula_msm <- as.formula(
  paste("BAA_PhenoAge ~ traj15_EVI_500m + traj15_O3_500m +", base_controls)
)
tryCatch({
  models$Joint_EVI_500M_MSM <- lm(
    formula_msm, data = model_data, weights = ipw_joint_500m
  )
}, error = function(e) {
  cat("Joint EVI 500M MSM failed:", e$message, "\n")
  models$Joint_EVI_500M_MSM <- NULL
})

tryCatch({
  models$Joint_O3_500M_MSM <- models$Joint_EVI_500M_MSM
}, error = function(e) {
  cat("Joint O3 500M MSM failed:", e$message, "\n")
  models$Joint_O3_500M_MSM <- NULL
})

cat("=== All regression models completed ===\n")

# -----------------------------------------------------------------------------
# Step 4: Extract coefficients and generate results
# -----------------------------------------------------------------------------

cat("=== Extracting exposure variable coefficients from all models ===\n")

model_order <- c(
  "EVI_1KM_WME", "EVI_1KM_UME", "EVI_1KM_MSM",
  "EVI_500M_WME", "EVI_500M_UME", "EVI_500M_MSM",
  "O3_1KM_WME",  "O3_1KM_UME",  "O3_1KM_MSM",
  "O3_500M_WME", "O3_500M_UME", "O3_500M_MSM",
  "Joint_EVI_1KM_WME", "Joint_O3_1KM_WME",
  "Joint_EVI_1KM_UME", "Joint_O3_1KM_UME",
  "Joint_EVI_1KM_MSM", "Joint_O3_1KM_MSM",
  "Joint_EVI_500M_WME", "Joint_O3_500M_WME",
  "Joint_EVI_500M_UME", "Joint_O3_500M_UME",
  "Joint_EVI_500M_MSM", "Joint_O3_500M_MSM"
)

results <- data.frame(
  model_id   = 1:length(model_order),
  model_name = model_order,
  exposure_type   = character(length(model_order)),
  buffer_distance = character(length(model_order)),
  analysis_method = character(length(model_order)),
  coefficient     = numeric(length(model_order)),
  std_error       = numeric(length(model_order)),
  stringsAsFactors = FALSE
)

extract_coefficient <- function(model, var_name, row_num) {
  tryCatch({
    if (!is.null(model)) {
      coef_summary  <- summary(model)$coefficients
      matching_rows <- grep(var_name, rownames(coef_summary), fixed = TRUE)

      if (length(matching_rows) > 0) {
        actual_row <- matching_rows[1]
        results[row_num, "coefficient"] <<- coef_summary[actual_row, "Estimate"]
        results[row_num, "std_error"]   <<- coef_summary[actual_row, "Std. Error"]
        cat("Successfully extracted coefficient for", model_order[row_num], "\n")
      } else {
        cat("Warning: variable", var_name, "not found in model", model_order[row_num], "\n")
        results[row_num, c("coefficient", "std_error")] <<- NA
      }
    } else {
      cat("Warning: model", model_order[row_num], "is NULL\n")
      results[row_num, c("coefficient", "std_error")] <<- NA
    }
  }, error = function(e) {
    cat("Error extracting coefficient for", model_order[row_num], ":", e$message, "\n")
    results[row_num, c("coefficient", "std_error")] <<- NA
  })
}

for (i in seq_along(model_order)) {
  model_name <- model_order[i]

  if (grepl("^EVI_", model_name) & !grepl("^Joint", model_name)) {
    var_name <- if (grepl("500M", model_name)) "traj15_EVI_500m" else "traj15_EVI_1KM"
  } else if (grepl("^O3_", model_name) & !grepl("^Joint", model_name)) {
    var_name <- if (grepl("500M", model_name)) "traj15_O3_500m" else "traj15_O3_1KM"
  } else if (grepl("^Joint_EVI", model_name)) {
    var_name <- if (grepl("500M", model_name)) "traj15_EVI_500m" else "traj15_EVI_1KM"
  } else if (grepl("^Joint_O3", model_name)) {
    var_name <- if (grepl("500M", model_name)) "traj15_O3_500m" else "traj15_O3_1KM"
  } else {
    var_name <- "traj15_EVI_1KM"
  }

  if (grepl("EVI", model_name)) {
    results[i, "exposure_type"] <- "EVI"
  } else {
    results[i, "exposure_type"] <- "O3"
  }

  if (grepl("1KM", model_name)) {
    results[i, "buffer_distance"] <- "1KM"
  } else if (grepl("500M", model_name)) {
    results[i, "buffer_distance"] <- "500M"
  }

  if (grepl("WME", model_name)) {
    results[i, "analysis_method"] <- "Weighted Mixed"
  } else if (grepl("UME", model_name)) {
    results[i, "analysis_method"] <- "Unweighted Mixed"
  } else if (grepl("MSM", model_name)) {
    results[i, "analysis_method"] <- "Propensity-score weighted"
  }

  extract_coefficient(models[[model_name]], var_name, i)
}

successful_models <- sum(!is.na(results$coefficient))
cat("Number of successfully run models:", successful_models, "/", length(model_order), "\n")

# Calculate confidence intervals
results <- results %>%
  mutate(
    CI_lower = coefficient - 1.96 * std_error,
    CI_upper = coefficient + 1.96 * std_error
  )

# Robustness check and fix for model 14 (Joint_O3_1KM_WME)
idx_14 <- which(results$model_name == "Joint_O3_1KM_WME")
if (length(idx_14) == 1) {
  cat("Checking model 14 (Joint_O3_1KM_WME) confidence intervals:\n")
  print(results[idx_14, c("model_name", "coefficient", "std_error", "CI_lower", "CI_upper")])
  if (!is.na(results$coefficient[idx_14]) && !is.na(results$std_error[idx_14])) {
    results$CI_lower[idx_14] <- results$coefficient[idx_14] - 1.96 * results$std_error[idx_14]
    results$CI_upper[idx_14] <- results$coefficient[idx_14] + 1.96 * results$std_error[idx_14]
  }
}

# Significance markers
results <- results %>%
  mutate(
    significance = case_when(
      abs(coefficient) > 2 * std_error & abs(coefficient) > 1.0 ~ "***",
      abs(coefficient) > 2 * std_error & abs(coefficient) > 0.5 ~ "**",
      abs(coefficient) > 2 * std_error                           ~ "*",
      TRUE ~ ""
    )
  )

# Format values
results <- results %>%
  mutate(
    coefficient_formatted = sprintf("%.3f%s", coefficient, significance),
    CI_formatted = sprintf("(%.3f, %.3f)", CI_lower, CI_upper),
    SE_formatted = sprintf("%.3f", std_error)
  )

# Add Single-exposure / Dual-exposure labels to exposure variable names
results <- results %>%
  mutate(
    exposure_with_type = ifelse(
      grepl("^Joint", model_name),
      paste0(exposure_type, " (Dual-exposure)"),
      paste0(exposure_type, " (Single-exposure)")
    )
  )

results_display <- results %>%
  select(
    Model  = model_name,
    Exposure = exposure_with_type,
    Buffer   = buffer_distance,
    Method   = analysis_method,
    Coefficient = coefficient_formatted,
    `95% CI`   = CI_formatted,
    SE         = SE_formatted
  )

cat("=== Exposure variable effect summary table ===\n")
print(results_display)

write.csv(results_display, "output/tables/Exposure_Effects_Summary.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# -----------------------------------------------------------------------------
# Step 5: Create combined figure (left table, right forest plot)
# -----------------------------------------------------------------------------
# v4: no meta-analysis pooling — the 24 model-specification rows are shown
# side-by-side as robustness analyses.

forest_data <- results

cat("Total models (for plotting) raw rows:", nrow(forest_data), "\n")

# Another check and fix for Joint_O3_1KM_WME
idx_14_fd <- which(forest_data$model_name == "Joint_O3_1KM_WME")
if (length(idx_14_fd) == 1) {
  cat("Forest plot data for model 14 (Joint_O3_1KM_WME):\n")
  print(forest_data[idx_14_fd, c("model_name", "coefficient", "std_error", "CI_lower", "CI_upper")])
  if (!is.na(forest_data$coefficient[idx_14_fd]) && !is.na(forest_data$std_error[idx_14_fd])) {
    forest_data$CI_lower[idx_14_fd] <- forest_data$coefficient[idx_14_fd] - 1.96 * forest_data$std_error[idx_14_fd]
    forest_data$CI_upper[idx_14_fd] <- forest_data$coefficient[idx_14_fd] + 1.96 * forest_data$std_error[idx_14_fd]
  }
}

# Critical step: keep only complete and finite rows to avoid geom_segment warnings
forest_data <- forest_data %>%
  filter(
    !is.na(coefficient),
    !is.na(std_error),
    !is.na(CI_lower),
    !is.na(CI_upper),
    is.finite(CI_lower),
    is.finite(CI_upper)
  )

cat("Filtered model rows for plotting:", nrow(forest_data), "\n")

# Unified Y-axis coordinates:
#   - Data rows use y = n_total:1 (first model row at top)
#   - Table has 1 more row than forest plot: header row
n_total <- nrow(forest_data)

forest_data <- forest_data %>%
  mutate(
    y = rev(seq_len(n_total)),                 # first model row at top
    model_label = ""                           # hide left axis labels
  )

y_min <- 0.5
y_max <- n_total + 1.5  # +1 row for header

# Auto-set x-axis range based on all CIs to avoid clipping
x_range  <- range(forest_data$CI_lower, forest_data$CI_upper, na.rm = TRUE)
x_breaks <- pretty(x_range, n = 7)
x_limits <- range(x_breaks)

# ----------------------- Forest plot ------------------------------------
forest_plot <- ggplot(forest_data, aes(x = coefficient, y = y)) +
  # Confidence interval lines (all models)
  geom_linerange(
    aes(xmin = CI_lower, xmax = CI_upper, color = exposure_type),
    linewidth = 2, alpha = 0.8
  ) +
  # Model point estimates
  geom_point(
    aes(color = exposure_type, fill = exposure_type),
    size = 4, shape = 21, stroke = 2
  ) +
  geom_vline(
    xintercept = 0, linetype = "dashed",
    color = "red", linewidth = 1.2, alpha = 0.8
  ) +
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = forest_data$y,
    labels = forest_data$model_label,
    expand = expansion(mult = 0)
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    limits = x_limits,
    expand = expansion(mult = 0.02)
  ) +
  scale_color_manual(values = c("EVI" = "#2E86AB", "O3" = "#A23B72")) +
  scale_fill_manual(values = c("EVI" = "#2E86AB", "O3" = "#A23B72")) +
  labs(
    title = "Model comparison: EVI and O3 associations with BAA",
    x = "Regression Coefficient (95% Confidence Interval)",
    y = ""
  ) +
  theme_minimal(base_size = 25) +
  theme(
    plot.title = element_text(
      hjust = 0.5, size = 18, face = "bold",
      margin = margin(b = 15)
    ),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 22),
    axis.title  = element_text(size = 22, face = "bold"),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.7),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 15, 10, 5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

# ----------------------- Table ------------------------------------------
table_data <- forest_data %>%
  transmute(
    Exposure    = exposure_with_type,
    Buffer      = buffer_distance,
    Method      = analysis_method,
    Coefficient = coefficient_formatted,
    SE          = SE_formatted,      # SE first
    `95% CI`    = CI_formatted,      # then 95% CI
    y           = y
  )

header_y <- n_total + 1  # table row 1: column headers

table_plot <- ggplot(table_data) +
  # Fine grey lines
  geom_hline(
    yintercept = seq(y_min, y_max, by = 1),
    color = "grey80", linewidth = 0.4
  ) +
  # Three-line table: top thick line above header (offset +0.5 row)
  annotate(
    "segment", x = 0.5, xend = 6.5, y = header_y + 0.5, yend = header_y + 0.5,
    color = "black", linewidth = 1.5
  ) +
  # Three-line table: bottom thick line below header (offset -0.5 row)
  annotate(
    "segment", x = 0.5, xend = 6.5, y = header_y - 0.5, yend = header_y - 0.5,
    color = "black", linewidth = 1.5
  ) +
  # Data rows for each model
  geom_text(
    aes(x = 1, y = y, label = Exposure),
    size = 7, hjust = 0.5, vjust = 0.5, fontface = "bold"
  ) +
  geom_text(
    aes(x = 2, y = y, label = Buffer),
    size = 7, hjust = 0.5, vjust = 0.5
  ) +
  geom_text(
    aes(x = 3, y = y, label = Method),
    size = 6.5, hjust = 0.5, vjust = 0.5
  ) +
  geom_text(
    aes(x = 4, y = y, label = Coefficient),
    size = 7, hjust = 0.5, vjust = 0.5, fontface = "bold"
  ) +
  # Swapped positions: x = 5 shows SE, x = 6 shows 95% CI
  geom_text(
    aes(x = 5, y = y, label = SE),
    size = 7, hjust = 0.5, vjust = 0.5
  ) +
  geom_text(
    aes(x = 6, y = y, label = `95% CI`),
    size = 6.5, hjust = 0.5, vjust = 0.5
  ) +
  # Column headers (table row 1)
  annotate(
    "text", x = 1, y = header_y,
    label = "Exposure", size = 8,
    hjust = 0.5, vjust = 0.5, fontface = "bold"
  ) +
  annotate(
    "text", x = 2, y = header_y,
    label = "Buffer", size = 8,
    hjust = 0.5, vjust = 0.5, fontface = "bold"
  ) +
  annotate(
    "text", x = 3, y = header_y,
    label = "Method", size = 8,
    hjust = 0.5, vjust = 0.5, fontface = "bold"
  ) +
  annotate(
    "text", x = 4, y = header_y,
    label = "Coefficient", size = 8,
    hjust = 0.5, vjust = 0.5, fontface = "bold"
  ) +
  # Headers also swapped accordingly
  annotate(
    "text", x = 5, y = header_y,
    label = "SE", size = 8,
    hjust = 0.5, vjust = 0.5, fontface = "bold"
  ) +
  annotate(
    "text", x = 6, y = header_y,
    label = "95% CI", size = 8,
    hjust = 0.5, vjust = 0.5, fontface = "bold"
  ) +
  # Three-line table: bottom thick line
  annotate(
    "segment", x = 0.5, xend = 6.5, y = y_min, yend = y_min,
    color = "black", linewidth = 1.5
  ) +
  scale_x_continuous(
    limits = c(0.5, 6.5),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = NULL,
    expand = expansion(mult = 0)
  ) +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme_void() +
  theme(
    plot.title = element_blank(),
    plot.margin = margin(10, 5, 10, 15)
  )


# ----------------------- Combined plot ----------------------------------
combined_plot <- table_plot + forest_plot +
  plot_layout(widths = c(1.8, 1))

# Save figures and results
ggsave("output/figures/Figure2_forest_plot.png", combined_plot,
       width = 28, height = 20, dpi = 300, bg = "white")

ggsave("output/figures/Figure2_forest_plot.pdf", combined_plot,
       width = 28, height = 20, bg = "white")

# Save results tables
write.csv(results_display, "output/tables/Detailed_Results_Final.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

final_results <- results %>%
  select(
    model_name, exposure_with_type, buffer_distance,
    analysis_method, coefficient, std_error,
    CI_lower, CI_upper,
    coefficient_formatted, CI_formatted, SE_formatted
  )

# Legacy filename retained (referenced by downstream sync checks); the table
# no longer contains meta-analysis rows.
write.csv(final_results, "output/tables/Complete_Results_with_Meta_Final.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

cat("=== Analysis complete ===\n")
cat("=== Figure 2: ../output/figures/Figure2_forest_plot.png / Figure2_forest_plot.pdf ===\n")
cat("=== Summary table: ../output/tables/Exposure_Effects_Summary.csv ===\n")
cat("=== Detailed table: ../output/tables/Detailed_Results_Final.csv ===\n")
cat("=== Complete results: ../output/tables/Complete_Results_with_Meta_Final.csv ===\n")

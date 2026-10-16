## ============================================================
## Quantile-based g-computation analysis (using qgcomp.glm.noboot)
## Exposures: EVI_1KM_mean4y and O3_1KM_mean4y
## Outputs:
##   1) Figure_qgcomp_quantile_gcomp_signed_weights.png
##      4-panel plot with beta_positive / beta_negative annotations at top
##   2) qgcomp_glm_boot_mixture_summary_EVI_O3_mean4y.csv
##      Columns: Exposure, Model Type, Buffer, Coefficient (SE), P-value, 95% CI
## ============================================================

library(dplyr)
library(ggplot2)
library(tibble)
library(haven)
library(qgcomp)

set.seed(20240101)

## -------------------------------
## 0. Working directory and data loading
## -------------------------------
data <- read_dta("DATA/analysis_ready/analysis_dataset.dta")

## ln(hh_consumption_total), if needed
if (!("ln_hh_consumption_total" %in% names(data)) && ("hh_consumption_total" %in% names(data))) {
  data <- data %>%
    mutate(ln_hh_consumption_total = ifelse(hh_consumption_total > 0, log(hh_consumption_total), NA_real_))
}

## -------------------------------
## 1. Climate variable PCA
## -------------------------------
climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")

climate_data <- data %>%
  select(all_of(climate_vars)) %>%
  na.omit()

pca_result <- prcomp(climate_data, scale. = TRUE)

climate_pcs <- as.data.frame(pca_result$x[, 1:5])
colnames(climate_pcs) <- paste0("climate_pc", 1:5)

## Simple column bind (consistent with earlier code)
data <- bind_cols(data, climate_pcs)

## -------------------------------
## 2. Exposures and covariates
## -------------------------------
mixture_vars <- c("EVI_1KM_mean4y", "O3_1KM_mean4y")

cov_age_sex <- intersect(c("age", "sex"), names(data))

cov_full <- intersect(
  c("age","sex","education","ln_hh_consumption_total","retirement","pension_insurance",
    "smoking_current","alcohol_current",
    "climate_pc1","climate_pc2","climate_pc3","climate_pc4","climate_pc5",
    "bmi","hypertension","dyslipidemia","year"),
  names(data)
)

vars_for_model <- unique(c("BAA_PhenoAge", mixture_vars, cov_full))

dat_qg <- data %>%
  select(any_of(vars_for_model)) %>%
  mutate(across(everything(), ~ ifelse(is.infinite(.), NA, .))) %>%
  na.omit()

cat("Sample size:", nrow(dat_qg), "\n")

## -------------------------------
## 3. Four model specifications (panels)
## -------------------------------
panel_specs <- list(
  list(
    name = "Unadjusted",
    covs = character(0),
    mode = "single"   # single-exposure model
  ),
  list(
    name = "Age- and sex-adjusted",
    covs = cov_age_sex,
    mode = "single"
  ),
  list(
    name = "Multivariable-adjusted (single)",
    covs = cov_full,
    mode = "single"
  ),
  list(
    name = "Multivariable-adjusted (dual)",
    covs = cov_full,
    mode = "dual"     # dual-exposure model
  )
)

## Note: Model Type for table output
mode_to_modeltype <- c(
  "single" = "Single-exposure",
  "dual"   = "Dual-exposure"
)

## Exposure labels
exposure_label_map <- c("EVI_1KM_mean4y" = "EVI",
                        "O3_1KM_mean4y"  = "O3")

panel_result    <- list()  # for plotting
mixture_summary <- list()  # for table output

## -------------------------------
## 4. Loop over panels
## -------------------------------
for (i in seq_along(panel_specs)) {
  spec  <- panel_specs[[i]]
  pname <- spec$name
  covs  <- spec$covs
  mode  <- spec$mode

  cat("=== Fitting panel:", pname, " (mode:", mode, ") ===\n")

  beta_vec  <- numeric(length(mixture_vars))
  se_vec    <- numeric(length(mixture_vars))
  p_vec     <- numeric(length(mixture_vars))
  ci_low    <- numeric(length(mixture_vars))
  ci_high   <- numeric(length(mixture_vars))
  names(beta_vec) <- names(se_vec) <- names(p_vec) <-
    names(ci_low) <- names(ci_high) <- mixture_vars

  if (mode == "single") {
    ## Single-exposure models: fit separate regressions for EVI and O3
    for (v in mixture_vars) {
      rhs <- if (length(covs) > 0) {
        paste(c(v, sprintf("`%s`", covs)), collapse = " + ")
      } else {
        v
      }
      f_qg <- as.formula(paste0("BAA_PhenoAge ~ ", rhs))

      fit <- qgcomp.glm.noboot(
        f_qg,
        data   = dat_qg,
        q      = 4,
        family = gaussian(),
        expnms = v
      )
      sm <- summary(fit$fit)
      coefs <- sm$coefficients

      est <- coefs[v, "Estimate"]
      se  <- coefs[v, "Std. Error"]
      p   <- coefs[v, ncol(coefs)]

      beta_vec[v] <- est
      se_vec[v]   <- se
      p_vec[v]    <- p
      ci_low[v]   <- est - 1.96 * se
      ci_high[v]  <- est + 1.96 * se
    }
  } else if (mode == "dual") {
    ## Dual-exposure model: both EVI and O3 in same regression
    rhs <- if (length(covs) > 0) {
      paste(c(mixture_vars, sprintf("`%s`", covs)), collapse = " + ")
    } else {
      paste(mixture_vars, collapse = " + ")
    }
    f_qg <- as.formula(paste0("BAA_PhenoAge ~ ", rhs))

    fit <- qgcomp.glm.noboot(
      f_qg,
      data   = dat_qg,
      q      = 4,
      family = gaussian(),
      expnms = mixture_vars
    )
    sm <- summary(fit$fit)
    coefs <- sm$coefficients

    for (v in mixture_vars) {
      est <- coefs[v, "Estimate"]
      se  <- coefs[v, "Std. Error"]
      p   <- coefs[v, ncol(coefs)]

      beta_vec[v] <- est
      se_vec[v]   <- se
      p_vec[v]    <- p
      ci_low[v]   <- est - 1.96 * se
      ci_high[v]  <- est + 1.96 * se
    }
  } else {
    stop("Unknown mode: ", mode)
  }

  ## 4.1 Signed weights (standardized beta)
  denom <- sum(abs(beta_vec))
  if (denom == 0) {
    w_signed <- rep(0, length(beta_vec))
    names(w_signed) <- names(beta_vec)
  } else {
    w_signed <- 100 * beta_vec / denom
  }

  ## 4.2 beta_positive and beta_negative (consistent with example figure)
  beta_positive <- sum(beta_vec[beta_vec > 0])
  beta_negative <- sum(abs(beta_vec[beta_vec < 0]))

  ## 4.3 Plotting data
  exposure_labels <- unname(exposure_label_map[mixture_vars])

  bars <- tibble(
    Panel    = pname,
    Exposure = factor(exposure_labels, levels = c("EVI","O3")),
    Signed_weight_percent = as.numeric(w_signed),
    label    = ifelse(
      Signed_weight_percent >= 0,
      paste0("+", sprintf("%.1f", abs(Signed_weight_percent)), "%"),
      paste0("-", sprintf("%.1f", abs(Signed_weight_percent)), "%")
    ),
    y_text   = ifelse(Signed_weight_percent >= 0,
                      Signed_weight_percent + 3.5,
                      Signed_weight_percent - 3.5)
  )

  y_top <- max(60, max(abs(bars$y_text)) + 8)

  ann <- tibble(
    Panel = pname,
    x     = 1.5,
    y     = y_top,
    label = sprintf("β_positive = %.3f\nβ_negative = %.3f",
                    beta_positive, beta_negative)
  )

  panel_result[[i]] <- list(bars = bars, ann = ann, ylim = y_top + 6)

  ## 4.4 Table data
  model_type <- mode_to_modeltype[mode]

  mixture_summary[[i]] <- tibble(
    Panel    = pname,
    ModelType = model_type,
    Exposure = exposure_labels,
    Buffer   = "1KM",
    Beta     = as.numeric(beta_vec),
    SE       = as.numeric(se_vec),
    P_value  = as.numeric(p_vec),
    Beta_CI_lower = as.numeric(ci_low),
    Beta_CI_upper = as.numeric(ci_high)
  )
}

## -------------------------------
## 5. Combine all panels for plotting
## -------------------------------
bars_all <- bind_rows(lapply(panel_result, function(x) x$bars))
ann_all  <- bind_rows(lapply(panel_result, function(x) x$ann))
y_lim    <- max(sapply(panel_result, function(x) x$ylim))

panel_order <- sapply(panel_specs, function(x) x$name)
bars_all$Panel <- factor(bars_all$Panel, levels = panel_order)
ann_all$Panel  <- factor(ann_all$Panel,  levels = panel_order)
bars_all$Exposure <- factor(bars_all$Exposure, levels = c("EVI","O3"))

legend_labels <- c("EVI", expression(O[3]))
pal <- c("EVI" = "#377EB8", "O3" = "#4DAF4A")

p_tripanel_vert <- ggplot(bars_all,
                          aes(x = Exposure,
                              y = Signed_weight_percent,
                              fill = Exposure)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             linewidth = 0.4, colour = "grey40") +
  geom_col(width = 0.7) +
  geom_text(aes(y = y_text, label = label),
            size = 3.8, fontface = "bold",
            vjust = ifelse(bars_all$Signed_weight_percent >= 0, 0, 1.1)) +
  geom_text(data = ann_all,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.4,
            vjust = 0, lineheight = 1.05) +
  facet_wrap(~ Panel, nrow = 1) +
  scale_fill_manual(values = pal,
                    breaks = c("EVI","O3"),
                    labels = legend_labels,
                    name   = NULL) +
  scale_x_discrete(labels = c("EVI", expression(O[3]))) +
  labs(
    x     = NULL,
    y     = "Signed weight (% of Σ|β|)",
    title = "Relative contribution of EVI and O₃ to the overall exposure composite effect\non biological age acceleration (quantile-based g-computation)"
  ) +
  coord_cartesian(ylim = c(-y_lim, y_lim), clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position    = "top",
    legend.justification = "left",
    plot.title         = element_text(face = "bold", hjust = 0),
    plot.margin        = unit(c(12, 24, 8, 8), "pt")
  )

print(p_tripanel_vert)

ggsave("output/figures/Figure_qgcomp_quantile_gcomp_signed_weights.png",
       p_tripanel_vert, width = 8.7, height = 4.4, dpi = 300)

## -------------------------------
## 6. Table output: same format as example
## -------------------------------
mixture_summary_df <- bind_rows(mixture_summary) %>%
  mutate(
    Beta          = round(Beta, 3),
    SE            = round(SE, 3),
    Beta_CI_lower = round(Beta_CI_lower, 3),
    Beta_CI_upper = round(Beta_CI_upper, 3),
    `Coefficient (SE)` = sprintf("%.3f (%.3f)", Beta, SE),
    P_text           = ifelse(P_value < 0.001, "<0.001",
                              sprintf("%.3f", P_value)),
    `95% CI`         = paste0("[",
                              sprintf("%.3f", Beta_CI_lower), ", ",
                              sprintf("%.3f", Beta_CI_upper), "]")
  )

## Output format:
## Exposure | Model Type | Buffer | Coefficient (SE) | P-value | 95% CI
table1_qgcomp <- mixture_summary_df %>%
  select(
    Exposure,
    `Model Type` = ModelType,
    Buffer,
    `Coefficient (SE)`,
    `P-value` = P_text,
    `95% CI`
  )

print(table1_qgcomp)

write.csv(table1_qgcomp,
          file = "output/tables/qgcomp_glm_boot_mixture_summary_EVI_O3_mean4y.csv",
          row.names = FALSE)

cat("\nOutput files:\n",
    "- ../output/figures/Figure_qgcomp_quantile_gcomp_signed_weights.png (4-panel, with beta_positive / beta_negative annotations)\n",
    "- ../output/tables/qgcomp_glm_boot_mixture_summary_EVI_O3_mean4y.csv (contains Single / Dual multivariable models)\n")

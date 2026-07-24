set.seed(20240101)

## ------------------------------------------------------------
## Quantile-based g-computation (qgcomp-style) for EVI / O3 mixture weights and regression coefficients
## Output for both 1KM and 500m buffer scales simultaneously
## ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(forcats)
library(tibble)
library(haven)

##-------------------------------------------------------------
## 0. Set working directory and read data
##-------------------------------------------------------------
data <- read_dta("DATA/analysis_ready/analysis_dataset.dta")

##-------------------------------------------------------------
## 1. Climate variable principal component analysis
##-------------------------------------------------------------
climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")

climate_data <- data %>%
  dplyr::select(all_of(climate_vars)) %>%
  na.omit()

pca_result <- prcomp(climate_data, scale. = TRUE)

# Extract first 5 principal components and merge into data
climate_pcs <- as.data.frame(pca_result$x[, 1:5])
colnames(climate_pcs) <- paste0("climate_pc", 1:5)
data <- bind_cols(data, climate_pcs)

##-------------------------------------------------------------
## 2. Quartile scoring function 0-3
##-------------------------------------------------------------
qscore <- function(x, q = 4){
  cut(
    x,
    breaks = quantile(x, probs = seq(0, 1, length = q + 1), na.rm = TRUE),
    include.lowest = TRUE,
    labels = FALSE
  ) - 1
}

##-------------------------------------------------------------
## 3. Wrapper function: run "qgcomp-style" analysis for a given buffer
##-------------------------------------------------------------
run_qgcomp_for_buffer <- function(data, evi_var, o3_var, buffer_label) {

  # Construct exposure quartile scores
  q_EVI_vec <- qscore(data[[evi_var]])
  q_O3_vec  <- qscore(data[[o3_var]])

  dat_q <- data %>%
    mutate(
      q_EVI = q_EVI_vec,
      q_O3  = q_O3_vec
    )

  # Covariate sets
  cov_age_sex <- intersect(c("age","sex"), names(dat_q))
  cov_full    <- intersect(
    c("age","sex","education","ln_hh_consumption_total","retirement","pension_insurance",
      "smoking_current","alcohol_current",
      "climate_pc1","climate_pc2","climate_pc3","climate_pc4","climate_pc5",
      "bmi","hypertension","dyslipidemia","year"),
    names(dat_q)
  )

  panels <- c("Unadjusted",
              "Age- and sex-adjusted",
              "Multivariable-adjusted (single)",
              "Multivariable-adjusted (dual)")

  panel_result <- list()
  model_list   <- list()   # store all models for this buffer

  ##------------------ First three panels: single-exposure models ------------------
  for(i in 1:3) {
    if(i == 1) {
      covs  <- character(0)
      pname <- "Unadjusted"
    } else if(i == 2) {
      covs  <- cov_age_sex
      pname <- "Age- and sex-adjusted"
    } else {
      covs  <- cov_full
      pname <- "Multivariable-adjusted (single)"
    }

    rhs <- if (length(covs) > 0) paste(sprintf("`%s`", covs), collapse = " + ") else "1"

    # Single-exposure: EVI
    f_lm_EVI <- as.formula(paste0("BAA_PhenoAge ~ q_EVI + ", rhs))
    mod_EVI  <- lm(f_lm_EVI, data = dat_q)

    # Single-exposure: O3
    f_lm_O3 <- as.formula(paste0("BAA_PhenoAge ~ q_O3 + ", rhs))
    mod_O3  <- lm(f_lm_O3,  data = dat_q)

    # Save models
    model_list[[paste0(pname, " - EVI")]] <- mod_EVI
    model_list[[paste0(pname, " - O3")]]  <- mod_O3

    # Extract only q_EVI / q_O3 coefficients for "signed weights"
    b_EVI <- coef(mod_EVI)["q_EVI"]
    b_O3  <- coef(mod_O3)["q_O3"]

    b <- c(EVI = b_EVI, O3 = b_O3)
    names(b) <- c("EVI","O3")

    sum_abs  <- sum(abs(b))
    w_signed <- 100 * b / sum_abs
    beta_pos <- sum(b[b > 0])
    beta_neg <- sum(abs(b[b < 0]))

    bars <- tibble(
      Buffer   = buffer_label,
      panel    = pname,
      exposure = factor(names(w_signed), levels = c("EVI","O3")),
      w_signed = as.numeric(w_signed),
      label    = ifelse(
        w_signed >= 0,
        paste0("+", sprintf("%.1f", abs(w_signed)), "%"),
        paste0("-", sprintf("%.1f", abs(w_signed)), "%")
      ),
      y_text   = ifelse(w_signed >= 0, w_signed + 3.5, w_signed - 3.5)
    )

    # Beta annotations placed 10 units above the tallest bar
    y_top <- max(bars$w_signed) + 10

    ann <- tibble(
      Buffer = buffer_label,
      panel  = pname,
      x      = 1.5,   # centered
      y      = y_top,
      label  = sprintf("beta_positive = %.3f\nbeta_negative = %.3f",
                       beta_pos, beta_neg)
    )

    panel_result[[i]] <- list(bars = bars, ann = ann)
  }

  ##------------------ Fourth panel: dual-exposure model ------------------
  covs  <- cov_full
  pname <- "Multivariable-adjusted (dual)"
  rhs   <- if (length(covs) > 0) paste(sprintf("`%s`", covs), collapse = " + ") else "1"

  f_lm <- as.formula(paste0("BAA_PhenoAge ~ q_EVI + q_O3 + ", rhs))
  mod  <- lm(f_lm, data = dat_q)

  model_list[[pname]] <- mod

  b_EVI <- coef(mod)["q_EVI"]
  b_O3  <- coef(mod)["q_O3"]

  b <- c(EVI = b_EVI, O3 = b_O3)
  names(b) <- c("EVI","O3")

  sum_abs  <- sum(abs(b))
  w_signed <- 100 * b / sum_abs
  beta_pos <- sum(b[b > 0])
  beta_neg <- sum(abs(b[b < 0]))

  bars <- tibble(
    Buffer   = buffer_label,
    panel    = pname,
    exposure = factor(names(w_signed), levels = c("EVI","O3")),
    w_signed = as.numeric(w_signed),
    label    = ifelse(
      w_signed >= 0,
      paste0("+", sprintf("%.1f", abs(w_signed)), "%"),
      paste0("-", sprintf("%.1f", abs(w_signed)), "%")
    ),
    y_text   = ifelse(w_signed >= 0, w_signed + 3.5, w_signed - 3.5)
  )

  y_top <- max(bars$w_signed) + 10

  ann <- tibble(
    Buffer = buffer_label,
    panel  = pname,
    x      = 1.5,
    y      = y_top,
    label  = sprintf("beta_positive = %.3f\nbeta_negative = %.3f",
                     beta_pos, beta_neg)
  )

  panel_result[[4]] <- list(bars = bars, ann = ann)

  ##------------------ Plotting data for this buffer ------------------
  bars_all <- bind_rows(lapply(panel_result, function(x) x$bars))
  ann_all  <- bind_rows(lapply(panel_result, function(x) x$ann))

  ##------------------ Extract coefficients / CIs / P values ------------------
  exposure_rows <- list()

  for(mn in names(model_list)) {
    m  <- model_list[[mn]]
    sm <- summary(m)
    coefs <- sm$coefficients
    ci    <- confint(m)

    is_dual <- grepl("dual", mn, ignore.case = TRUE)

    if (is_dual) {
      terms <- intersect(c("q_EVI", "q_O3"), rownames(coefs))
    } else {
      if (grepl("EVI", mn)) {
        terms <- intersect("q_EVI", rownames(coefs))
      } else if (grepl("O3", mn)) {
        terms <- intersect("q_O3", rownames(coefs))
      } else {
        next
      }
    }

    for(term in terms) {
      row <- coefs[term, ]
      est <- row["Estimate"]
      se  <- row["Std. Error"]
      t_v <- row["t value"]
      p_v <- row["Pr(>|t|)"]
      ci_l <- ci[term, 1]
      ci_u <- ci[term, 2]

      exposure_name <- ifelse(term == "q_EVI", "EVI", "O3")
      model_type    <- ifelse(is_dual, "Dual-exposure", "Single-exposure")

      exposure_rows[[length(exposure_rows) + 1]] <- data.frame(
        Buffer     = buffer_label,
        Exposure   = exposure_name,
        ModelType  = model_type,
        ModelLabel = mn,
        Estimate   = est,
        StdError   = se,
        t_value    = t_v,
        P_value    = p_v,
        CI_lower   = ci_l,
        CI_upper   = ci_u,
        stringsAsFactors = FALSE
      )
    }
  }

  exposure_summary <- bind_rows(exposure_rows) %>%
    mutate(
      `Coefficient (SE)` = sprintf("%.3f (%.3f)", Estimate, StdError),
      P_text             = ifelse(P_value < 0.001, "<0.001",
                                  sprintf("%.3f", P_value)),
      `95% CI`           = paste0("[",
                                  sprintf("%.3f", CI_lower), ", ",
                                  sprintf("%.3f", CI_upper), "]")
    )

  table_to_export <- exposure_summary %>%
    dplyr::select(
      Buffer,
      Exposure,
      `Model Type` = ModelType,
      `Coefficient (SE)`,
      `P-value` = P_text,
      `95% CI`
    )

  list(
    bars_all   = bars_all,
    ann_all    = ann_all,
    coef_table = table_to_export
  )
}

##-------------------------------------------------------------
## 4. Run "qgcomp-style" analysis for 1KM and 500m buffers separately
##-------------------------------------------------------------
res_1KM  <- run_qgcomp_for_buffer(
  data,
  evi_var = "EVI_1KM_mean4y",
  o3_var  = "O3_1KM_mean4y",
  buffer_label = "1KM"
)

res_500m <- run_qgcomp_for_buffer(
  data,
  evi_var = "EVI_500m_mean4y",
  o3_var  = "O3_500m_mean4y",
  buffer_label = "500m"
)

##-------------------------------------------------------------
## 5. Merge plotting data from both scales and plot
##-------------------------------------------------------------
bars_all_combined <- bind_rows(res_1KM$bars_all, res_500m$bars_all)
ann_all_combined  <- bind_rows(res_1KM$ann_all,  res_500m$ann_all)

panel_order <- c("Unadjusted",
                 "Age- and sex-adjusted",
                 "Multivariable-adjusted (single)",
                 "Multivariable-adjusted (dual)")

bars_all_combined$panel    <- factor(bars_all_combined$panel, levels = panel_order)
ann_all_combined$panel     <- factor(ann_all_combined$panel,  levels = panel_order)
bars_all_combined$Buffer   <- factor(bars_all_combined$Buffer, levels = c("1KM","500m"))
ann_all_combined$Buffer    <- factor(ann_all_combined$Buffer,  levels = c("1KM","500m"))
bars_all_combined$exposure <- factor(bars_all_combined$exposure, levels = c("EVI","O3"))

## -- Key: separately compute upper and lower y-axis ranges --
top_lim <- max(c(
  bars_all_combined$w_signed[bars_all_combined$w_signed > 0],
  bars_all_combined$y_text[bars_all_combined$y_text > 0],
  ann_all_combined$y
)) + 8   # extra top space to avoid clipping 500m beta annotations

bottom_lim <- max(abs(c(
  bars_all_combined$w_signed[bars_all_combined$w_signed < 0],
  bars_all_combined$y_text[bars_all_combined$y_text < 0]
))) + 5   # bottom space kept similar to original

legend_labels <- c("EVI", expression(O[3]))
pal <- c("EVI" = "#377EB8", "O3" = "#4DAF4A")

panel_labs <- c(
  "Unadjusted"                      = "Unadjusted",
  "Age- and sex-adjusted"           = "Age- and\nsex-adjusted",
  "Multivariable-adjusted (single)" = "Multivariable-adjusted\n(single)",
  "Multivariable-adjusted (dual)"   = "Multivariable-adjusted\n(dual)"
)

p_combined <- ggplot(bars_all_combined,
                     aes(x = exposure, y = w_signed, fill = exposure)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  geom_col(width = 0.7) +
  geom_text(aes(y = y_text, label = label),
            size = 3.8, fontface = "bold",
            vjust = ifelse(bars_all_combined$w_signed >= 0, 0, 1.1)) +
  geom_text(data = ann_all_combined,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.6,
            vjust = 0, lineheight = 1.05) +
  facet_grid(Buffer ~ panel,
             labeller = labeller(panel = panel_labs)) +
  scale_fill_manual(values = pal,
                    breaks = c("EVI","O3"),
                    labels = legend_labels,
                    name   = NULL) +
  labs(
    x     = NULL,
    y     = "Signed weight (% of Sigma|beta|)",
    title = "EVI and O3: Signed weights of joint exposure by buffer scale based on quantile g-computation"
  ) +
  # Use asymmetric y-axis range: top raised, bottom mostly unchanged
  coord_cartesian(ylim = c(-bottom_lim, top_lim)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x  = element_blank(),
    legend.position     = "top",
    legend.justification = "left",
    plot.title          = element_text(face = "bold"),
    plot.margin         = unit(c(12, 24, 8, 8), "pt"),
    strip.text.x        = element_text(size = 9,
                                       margin = margin(t = 3, r = 0, b = 3, l = 0)),
    panel.spacing.y     = unit(5, "pt")   # keep spacing between upper and lower panel rows
  )

print(p_combined)

ggsave("output/figures/qgcomp_1KM_500m_signed_weights.png",
       p_combined, width = 6.7, height = 8, dpi = 300)


##-------------------------------------------------------------
## 6. Merge regression coefficient tables from both scales and output
##-------------------------------------------------------------
coef_all_buffers <- bind_rows(res_1KM$coef_table, res_500m$coef_table)

print(coef_all_buffers)

write.csv(coef_all_buffers,
          file = "output/tables/gcomp_core_exposures_summary_1KM_500m.csv",
          row.names = FALSE)

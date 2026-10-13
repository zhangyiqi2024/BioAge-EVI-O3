set.seed(20240101)

library(dplyr)
library(ggplot2)
library(forcats)
library(tibble)
library(haven)

#-------------------------------
# 0. Set working directory and read data
#-------------------------------
data <- read_dta("DATA/analysis_ready/analysis_dataset.dta")

#-------------------------------
# 1. Climate variable PCA
#-------------------------------
climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")
climate_data <- data %>%
  select(all_of(climate_vars)) %>%
  na.omit()

pca_result <- prcomp(climate_data, scale. = TRUE)

# Extract first 5 principal components
climate_pcs <- as.data.frame(pca_result$x[, 1:5])
colnames(climate_pcs) <- paste0("climate_pc", 1:5)

# Merge back to original data
data <- bind_cols(data, climate_pcs)

#-------------------------------
# 2. Quartile scoring 0-3
#-------------------------------
qscore <- function(x, q = 4){
  cut(
    x,
    breaks = quantile(x, probs = seq(0, 1, length = q + 1), na.rm = TRUE),
    include.lowest = TRUE,
    labels = FALSE
  ) - 1
}

dat_q <- data %>%
  mutate(
    q_EVI = qscore(`EVI_1KM_mean4y`),
    q_O3  = qscore(`O3_1KM_mean4y`)
  )

#-------------------------------
# 3. Covariate sets
#-------------------------------
cov_age_sex <- intersect(c("age","sex"), names(dat_q))
cov_full    <- intersect(
  c("age","sex","education","ln_hh_consumption_total","retirement","pension_insurance",
    "smoking_current","alcohol_current",
    "climate_pc1","climate_pc2","climate_pc3","climate_pc4","climate_pc5",
    "bmi","hypertension","dyslipidemia","year"),
  names(dat_q)
)

#-------------------------------
# 4. Fit models and prepare plotting data
#-------------------------------
panels <- c("Unadjusted",
            "Age- and sex-adjusted",
            "Multivariable-adjusted (single)",
            "Multivariable-adjusted (dual)")

panel_result <- list()

# Save all regression models for later coefficient/CI/P-value extraction
model_list <- list()

#---- First three panels: single-exposure models ----
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

  # Extract only q_EVI / q_O3 coefficients for signed weights
  b_EVI <- coef(mod_EVI)["q_EVI"]
  b_O3  <- coef(mod_O3)["q_O3"]

  b <- c(EVI = b_EVI, O3 = b_O3)
  names(b) <- c("EVI","O3")

  sum_abs  <- sum(abs(b))
  w_signed <- 100 * b / sum_abs
  beta_pos <- sum(b[b > 0])
  beta_neg <- sum(abs(b[b < 0]))

  bars <- tibble(
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

  y_top <- max(60, max(abs(bars$y_text)) + 8)
  ann <- tibble(
    panel = pname,
    x     = 1.5,
    y     = y_top,
    label = sprintf("beta_positive = %.3f\nbeta_negative = %.3f",
                    beta_pos, beta_neg)
  )

  panel_result[[i]] <- list(bars = bars, ann = ann, ylim = y_top + 6)
}

#---- Fourth panel: dual-exposure model ----
covs  <- cov_full
pname <- "Multivariable-adjusted (dual)"
rhs   <- if (length(covs) > 0) paste(sprintf("`%s`", covs), collapse = " + ") else "1"

f_lm <- as.formula(paste0("BAA_PhenoAge ~ q_EVI + q_O3 + ", rhs))
mod  <- lm(f_lm, data = dat_q)

# Save dual-exposure model
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

y_top <- max(60, max(abs(bars$y_text)) + 8)
ann <- tibble(
  panel = pname,
  x     = 1.5,
  y     = y_top,
  label = sprintf("beta_positive = %.3f\nbeta_negative = %.3f",
                  beta_pos, beta_neg)
)

panel_result[[4]] <- list(bars = bars, ann = ann, ylim = y_top + 6)

#-------------------------------
# 5. Combine all panels for plotting
#-------------------------------
bars_all <- bind_rows(lapply(panel_result, function(x) x$bars))
ann_all  <- bind_rows(lapply(panel_result, function(x) x$ann))
y_lim    <- max(sapply(panel_result, function(x) x$ylim))

panel_order <- c("Unadjusted",
                 "Age- and sex-adjusted",
                 "Multivariable-adjusted (single)",
                 "Multivariable-adjusted (dual)")

bars_all$panel   <- factor(bars_all$panel, levels = panel_order)
ann_all$panel    <- factor(ann_all$panel,  levels = panel_order)
bars_all$exposure <- factor(bars_all$exposure, levels = c("EVI","O3"))

# Legend and color palette
legend_labels <- c("EVI", expression(O[3]))
pal <- c("EVI" = "#377EB8", "O3" = "#4DAF4A")

#-------------------------------
# 6. Plot
#-------------------------------
p_tripanel_vert <- ggplot(bars_all, aes(x = exposure, y = w_signed, fill = exposure)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  geom_col(width = 0.7) +
  geom_text(aes(y = y_text, label = label),
            size = 3.8, fontface = "bold",
            vjust = ifelse(bars_all$w_signed >= 0, 0, 1.1)) +
  geom_text(data = ann_all, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.6, vjust = 0, lineheight = 1.05) +
  facet_wrap(~panel, nrow = 1) +
  scale_fill_manual(values = pal,
                    breaks = c("EVI","O3"),
                    labels = legend_labels,
                    name   = NULL) +
  # X-axis tick labels: EVI and O3 (3 as subscript)
  scale_x_discrete(labels = c("EVI", expression(O[3]))) +
  labs(
    x     = NULL,
    y     = "Signed weight (% of Sigma|beta|)",
    title = "EVI and O3: Signed weights of joint exposure"
  ) +
  coord_cartesian(ylim = c(-y_lim, y_lim), clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position    = "top",
    legend.justification = "left",
    plot.title         = element_text(face = "bold"),
    plot.margin        = unit(c(12, 24, 8, 8), "pt")
  )


print(p_tripanel_vert)

# Canvas width reduced by 1/3: 13 * 2/3 approx 8.7
ggsave("output/figures/qgcomp_signed_weights.png", p_tripanel_vert,
       width = 8.7, height = 4.4, dpi = 300)

#-------------------------------
# 7. Core exposure variable coefficients, 95% CI and P value summary table
#-------------------------------
# Extract only q_EVI and q_O3 coefficients (core exposure variables) and output to working directory

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
    buffer_name   <- "1KM"

    exposure_rows[[length(exposure_rows) + 1]] <- data.frame(
      Exposure   = exposure_name,
      ModelType  = model_type,
      Buffer     = buffer_name,
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

# Full numeric version for supplementary materials / further analysis
print(exposure_summary)

# Export a journal-style concise table
table_to_export <- exposure_summary %>%
  select(
    Exposure,
    `Model Type` = ModelType,
    Buffer,
    `Coefficient (SE)`,
    `P-value` = P_text,
    `95% CI`
  )

write.csv(table_to_export,
          file = "output/tables/gcomp_core_exposures_summary.csv",
          row.names = FALSE)

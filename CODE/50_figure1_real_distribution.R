# =============================================================================
# 50_figure1_real_distribution.R  (v5 revision, task T11a)
# New main Figure 1: OBSERVED BAA distributions across the four joint
# EVI x O3 exposure groups (violin + boxplot), 1-km and 500-m buffers.
# Replaces the model-derived boxplots of v4 with real data distributions.
# White points with error bars overlay group means +/- 95% CI of the mean
# (corrective pass 2026-07-24, user request).
# Okabe-Ito colourblind-safe palette, Nature double-column specs.
# Input: DATA/analysis_ready/analysis_dataset.dta
# Also writes the figure source data (03_source_data schema, with CI columns)
# to output/tables/figure1_observed_distribution_v5.csv
# =============================================================================

suppressPackageStartupMessages({
  library(haven); library(dplyr); library(ggplot2); library(patchwork)
})

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
sink("logs/50_figure1_real_distribution.log", split = TRUE)
on.exit(sink(), add = TRUE)

theme_nature <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line      = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks     = element_line(linewidth = 0.35, colour = "black"),
      axis.title     = element_text(size = base_size),
      axis.text      = element_text(size = base_size - 0.5, colour = "#272727"),
      strip.background = element_blank(),
      strip.text     = element_text(size = base_size, face = "bold", hjust = 0),
      plot.tag       = element_text(size = base_size + 2, face = "bold"),
      panel.grid     = element_blank(),
      legend.position = "none"
    )
}
theme_set(theme_nature())

group_palette <- c(
  "Low EVI / low O3"   = "#0072B2",
  "Low EVI / high O3"  = "#E69F00",
  "High EVI / low O3"  = "#009E73",
  "High EVI / high O3" = "#CC79A7")
group_levels <- names(group_palette)
x_labels <- c(
  "Low EVI / low O3"   = expression(atop("Low EVI", "low O"[3])),
  "Low EVI / high O3"  = expression(atop("Low EVI", "high O"[3])),
  "High EVI / low O3"  = expression(atop("High EVI", "low O"[3])),
  "High EVI / high O3" = expression(atop("High EVI", "high O"[3])))

df <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
df <- df[df$age > 45 & df$age < 90, ]

lab_map <- c("1" = "Low EVI / low O3", "2" = "Low EVI / high O3",
             "3" = "High EVI / low O3", "4" = "High EVI / high O3")

# per-group summary statistics (observed): mean, sd, quartiles and the 95% CI
# of the mean (mean +/- 1.96 * sd / sqrt(n)) for the overlay layer
group_stats <- function(d, gvar, scale_lab) {
  d <- d[!is.na(d[[gvar]]) & !is.na(d$BAA_PhenoAge), ]
  d$group <- factor(lab_map[as.character(as.numeric(d[[gvar]]))],
                    levels = group_levels)
  d <- d[!is.na(d$group), ]
  d %>% group_by(group) %>% summarise(
    n = n(), mean = mean(BAA_PhenoAge), sd = sd(BAA_PhenoAge),
    median = median(BAA_PhenoAge),
    q1 = quantile(BAA_PhenoAge, 0.25, names = FALSE),
    q3 = quantile(BAA_PhenoAge, 0.75, names = FALSE),
    .groups = "drop") %>%
    mutate(scale = scale_lab, se = sd / sqrt(n),
           ci_lower = mean - 1.96 * se, ci_upper = mean + 1.96 * se) %>%
    select(scale, group, n, mean, sd, median, q1, q3, se, ci_lower, ci_upper)
}

stats_all <- bind_rows(group_stats(df, "joint_exp_1KM", "1KM"),
                       group_stats(df, "joint_exp_500m", "500m"))
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(as.data.frame(stats_all),
          "output/tables/figure1_observed_distribution_v5.csv", row.names = FALSE)
print(as.data.frame(stats_all), row.names = FALSE)

make_panel <- function(d, gvar, tag) {
  d <- d[!is.na(d[[gvar]]) & !is.na(d$BAA_PhenoAge), ]
  d$group <- factor(lab_map[as.character(as.numeric(d[[gvar]]))],
                    levels = group_levels)
  d <- d[!is.na(d$group), ]
  ns <- d %>% count(group)
  st <- stats_all %>% filter(scale == ifelse(gvar == "joint_exp_1KM", "1KM", "500m"))
  cat(tag, "N per group:", paste(ns$group, ns$n, sep = "=", collapse = ", "), "\n")
  ggplot(d, aes(x = group, y = BAA_PhenoAge, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.55, colour = NA, scale = "width") +
    geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9,
                 colour = "black", linewidth = 0.35) +
    # mean +/- 95% CI of the mean: white-filled point with black border + error bar
    geom_errorbar(data = st, aes(x = group, ymin = ci_lower, ymax = ci_upper),
                  inherit.aes = FALSE, width = 0.10, linewidth = 0.45,
                  colour = "black") +
    geom_point(data = st, aes(x = group, y = mean), inherit.aes = FALSE,
               shape = 21, size = 2.2, fill = "white", colour = "black",
               stroke = 0.6) +
    geom_text(data = ns, aes(x = group, y = Inf, label = paste0("n = ", n)),
              vjust = 1.6, size = 2.2, colour = "#272727", inherit.aes = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "#D8D8D8", linewidth = 0.35) +
    scale_fill_manual(values = group_palette) +
    scale_x_discrete(labels = x_labels) +
    labs(x = NULL, y = "Observed BAA (years)", tag = tag) +
    theme(axis.ticks.x = element_blank())
}

p_a <- make_panel(df, "joint_exp_1KM", "a") +
  labs(subtitle = "1 km residential buffer") +
  theme(plot.subtitle = element_text(size = 7, hjust = 0))
p_b <- make_panel(df, "joint_exp_500m", "b") +
  labs(subtitle = "500 m residential buffer") +
  theme(plot.subtitle = element_text(size = 7, hjust = 0))

fig1 <- p_a / p_b

save_all <- function(plot, stem, width_mm = 183, height_mm = 190) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h,
                       family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = 600, compression = "lzw"); print(plot); dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = w, height = h, units = "in",
                res = 300); print(plot); dev.off()
}
save_all(fig1, "output/figures/Figure1_observed_distribution")
cat("\nDone. output/figures/Figure1_observed_distribution.{svg,pdf,tiff,png}\n")

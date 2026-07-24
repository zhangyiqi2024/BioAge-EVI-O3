# =============================================================================
# 51_figure2_continuous_forest.R  (v5 revision, task T11b)
# New main Figure 2: PRIMARY continuous dual-exposure associations.
# Per-IQR adjusted differences in BAA for EVI and O3, 1-km and 500-m buffers,
# unweighted (primary) and IPW (sensitivity), from WeMix mixed models with
# robust variance (43_primary_continuous_models.R).
# Okabe-Ito palette; Nature double-column specs.
# Input: output/tables/primary_continuous_models_v5.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork)
})

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
sink("logs/51_figure2_continuous_forest.log", split = TRUE)
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
      panel.grid     = element_blank(),
      legend.position = "bottom",
      legend.title   = element_text(size = base_size - 0.5),
      legend.text    = element_text(size = base_size - 0.8)
    )
}
theme_set(theme_nature())

d <- read.csv("output/tables/primary_continuous_models_v5.csv",
              stringsAsFactors = FALSE) %>%
  filter(model == "dual-exposure") %>%
  mutate(
    exposure = factor(exposure, levels = c("EVI", "O3")),
    scale    = factor(scale, levels = c("1KM", "500m")),
    spec     = factor(weight_spec, levels = c("unweighted", "ipw"),
                      labels = c("Unweighted (primary)", "IPW (sensitivity)")),
    anno     = sprintf("%+.2f (%+.2f, %+.2f), P = %.3f",
                       beta_iqr, lci_iqr, uci_iqr, p),
    row_lab  = interaction(scale, spec, sep = " — "))

cat("Rows:", nrow(d), "\n"); print(d[, c("exposure", "scale", "spec",
      "beta_iqr", "lci_iqr", "uci_iqr", "p")])

spec_col <- c("Unweighted (primary)" = "#0072B2", "IPW (sensitivity)" = "#E69F00")
spec_pch <- c("Unweighted (primary)" = 16, "IPW (sensitivity)" = 17)

make_panel <- function(dd, ylab, tag) {
  dd <- dd %>% arrange(scale, desc(spec)) %>%
    mutate(y = row_number())
  ggplot(dd, aes(x = beta_iqr, y = y, colour = spec, shape = spec)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "#767676", linewidth = 0.35) +
    geom_errorbar(aes(xmin = lci_iqr, xmax = uci_iqr),
                  orientation = "y", width = 0.22, linewidth = 0.55) +
    geom_point(size = 2.4, stroke = 0.35) +
    geom_text(aes(x = uci_iqr + 0.12, label = anno),
              hjust = 0, size = 2.1, colour = "#272727", show.legend = FALSE) +
    scale_y_continuous(breaks = dd$y, labels = dd$row_lab,
                       trans = "reverse", expand = expansion(add = 0.6)) +
    scale_colour_manual(values = spec_col, name = NULL) +
    scale_shape_manual(values = spec_pch, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.42))) +
    labs(x = ylab, y = NULL, tag = tag) +
    theme(panel.grid.major.x = element_line(colour = "#ECECEC",
                                            linewidth = 0.3))
}

p_evi <- make_panel(d %>% filter(exposure == "EVI"),
                    "Adjusted BAA difference per IQR higher EVI (years)", "a")
p_o3  <- make_panel(d %>% filter(exposure == "O3"),
                    expression("Adjusted BAA difference per IQR higher O"[3]*" (years)"),
                    "b")

fig2 <- p_evi / p_o3 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_all <- function(plot, stem, width_mm = 183, height_mm = 120) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h,
                       family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = 600, compression = "lzw"); print(plot); dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = w, height = h, units = "in",
                res = 300); print(plot); dev.off()
}
save_all(fig2, "output/figures/Figure2_continuous_forest")
cat("\nDone. output/figures/Figure2_continuous_forest.{svg,pdf,tiff,png}\n")

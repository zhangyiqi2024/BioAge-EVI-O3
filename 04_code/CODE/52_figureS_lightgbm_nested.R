# =============================================================================
# 52_figureS_lightgbm_nested.R  (v5 revision, task T11c)
# New SI figure: participant-level split LightGBM re-evaluation (49_lightgbm_nested.R).
#   a: grouped exact TreeSHAP shares on the held-out test set
#   b: permutation importance (test-set R2 drop) for the top features
# Annotated with the honest held-out test R2.
# Inputs: output/tables/lightgbm_nested_{performance,shap_groups,importance}_v5.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork)
})

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
sink("logs/52_figureS_lightgbm_nested.log", split = TRUE)
on.exit(sink(), add = TRUE)

theme_nature <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line      = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks     = element_line(linewidth = 0.35, colour = "black"),
      axis.title     = element_text(size = base_size),
      axis.text      = element_text(size = base_size - 0.5, colour = "#272727"),
      strip.background = element_blank(),
      plot.tag       = element_text(size = base_size + 2, face = "bold"),
      panel.grid     = element_blank(),
      legend.position = "none"
    )
}
theme_set(theme_nature())

perf <- read.csv("output/tables/lightgbm_nested_performance_v5.csv")
gsh  <- read.csv("output/tables/lightgbm_nested_shap_groups_v5.csv")
imp  <- read.csv("output/tables/lightgbm_nested_importance_v5.csv")

anno <- sprintf("Held-out test R² = %.3f (participant-level split)", perf$r2_test)
cat(anno, "\n")

p_a <- ggplot(gsh, aes(x = reorder(group, share_pct), y = share_pct)) +
  geom_col(fill = "#0072B2", alpha = 0.85, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", share_pct)),
            hjust = -0.1, size = 2.2, colour = "#272727") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Share of mean |SHAP| (%)", tag = "a",
       caption = anno) +
  theme(plot.caption = element_text(size = 6.5, hjust = 0, colour = "#767676"))

top <- imp %>% arrange(desc(r2_drop_mean)) %>% head(10)
p_b <- ggplot(top, aes(x = reorder(feature, r2_drop_mean), y = r2_drop_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "#D8D8D8", linewidth = 0.35) +
  geom_col(fill = "#E69F00", alpha = 0.85, width = 0.7) +
  coord_flip() +
  labs(x = NULL, y = "Permutation importance (test-set R² drop)", tag = "b")

figS <- p_a | p_b

save_all <- function(plot, stem, width_mm = 183, height_mm = 90) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h,
                       family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = 600, compression = "lzw"); print(plot); dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = w, height = h, units = "in",
                res = 300); print(plot); dev.off()
}
save_all(figS, "output/figures/FigS8_lightgbm_nested")
cat("\nDone. output/figures/FigS8_lightgbm_nested.{svg,pdf,tiff,png}\n")

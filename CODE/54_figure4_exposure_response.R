# =============================================================================
# 54_figure4_exposure_response.R  (v5 corrective pass, 2026-07-24, user request)
# New main Figure 4: exposure-response associations of EVI and O3 with BAA.
#   a, b: average predicted BAA across the observed range of each 4-year mean
#         exposure (natural cubic splines, df = 3) under marginal
#         standardization over the primary analytical frame; 95% CI ribbons;
#         primary linear specification overlaid (dashed); IQR band + rug.
#   c:    average predicted BAA over the joint EVI x O3 plane (linear dual
#         model, consistent with the primary specification); observed
#         participants overlaid; grid cells far from the data cloud masked.
# All predictions use the robust (sandwich/Taylor) covariance from WeMix
# (same estimator as the primary analysis, script 43). Descriptive
# association figure; no causal or offsetting claims.
# Input:  DATA/analysis_ready/analysis_dataset.dta
# Output: output/figures/Figure4_exposure_response.{svg,pdf,tiff,png}
#         output/tables/figure4_exposure_response_v5.csv (curve source data)
#         output/tables/figure4_joint_contour_v5.csv     (contour source data)
# =============================================================================

suppressPackageStartupMessages({
  library(haven); library(dplyr); library(ggplot2); library(patchwork)
  library(splines); library(WeMix)
})

set.seed(20240101)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
sink("logs/54_figure4_exposure_response.log", split = TRUE)
on.exit(sink(), add = TRUE)

theme_nature <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line  = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text  = element_text(size = base_size - 0.5, colour = "#272727"),
      plot.tag   = element_text(size = base_size + 2, face = "bold"),
      panel.grid = element_blank())
}

# ---- data: primary analytical frame (script 43 convention) -------------------
df <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
df <- df[df$age > 45 & df$age < 90, ]
df$ln_hh_consumption_total <- ifelse(df$hh_consumption_total > 0,
                                     log(df$hh_consumption_total), NA_real_)
climate_vars <- c("TMP", "PER", "SNOW", "DEM", "lat")
cc <- complete.cases(df[, climate_vars])
pc_fit <- prcomp(df[cc, climate_vars], center = TRUE, scale. = TRUE)
df[, paste0("climate_pc", 1:5)] <- predict(pc_fit, newdata = df[, climate_vars])[, 1:5]

covars <- c("age", "sex", "education", "ln_hh_consumption_total",
            "retirement", "pension_insurance", "smoking_current",
            "alcohol_current", paste0("climate_pc", 1:5),
            "bmi", "hypertension", "dyslipidemia", "year")
mvars <- c("BAA_PhenoAge", "EVI_1KM_mean4y", "O3_1KM_mean4y", "id_num", covars)
md <- df[complete.cases(df[, mvars]), mvars]
md$id_num <- factor(md$id_num); md$year <- factor(md$year)
for (v in c("sex", "education", "retirement", "pension_insurance",
            "smoking_current", "alcohol_current", "hypertension", "dyslipidemia"))
  md[[v]] <- factor(md[[v]])
md$EVIx <- md$EVI_1KM_mean4y; md$O3x <- md$O3_1KM_mean4y
cat("Primary frame:", nrow(md), "obs,", length(unique(md$id_num)), "participants\n")

# ---- spline bases (knots stored for out-of-grid prediction) ------------------
df_spline <- 3
B_e <- ns(md$EVIx, df = df_spline)
B_o <- ns(md$O3x,  df = df_spline)
colnames(B_e) <- paste0("nsEVI", 1:df_spline)
colnames(B_o) <- paste0("nsO3",  1:df_spline)
sdat <- cbind(md, B_e, B_o)
sdat$w <- 1; sdat$w2 <- 1

f_ns <- as.formula(paste(
  "BAA_PhenoAge ~ nsEVI1 + nsEVI2 + nsEVI3 + nsO31 + nsO32 + nsO33 +",
  "age + sex + education + ln_hh_consumption_total + retirement +",
  "pension_insurance + smoking_current + alcohol_current +",
  "climate_pc1 + climate_pc2 + climate_pc3 + climate_pc4 + climate_pc5 +",
  "bmi + hypertension + dyslipidemia + year"))
fit_ns <- mix(update(f_ns, . ~ . + (1 | id_num)), data = sdat,
              weights = c("w", "w2"), cWeights = FALSE, verbose = FALSE)
beta_ns <- coef(fit_ns); V_ns <- fit_ns$cov_mat; beta_ns <- beta_ns[rownames(V_ns)]

f_lin <- BAA_PhenoAge ~ EVIx + O3x + age + sex + education +
  ln_hh_consumption_total + retirement + pension_insurance +
  smoking_current + alcohol_current + climate_pc1 + climate_pc2 +
  climate_pc3 + climate_pc4 + climate_pc5 + bmi + hypertension +
  dyslipidemia + year
fit_lin <- mix(update(f_lin, . ~ . + (1 | id_num)), data = sdat,
               weights = c("w", "w2"), cWeights = FALSE, verbose = FALSE)
beta_ln <- coef(fit_lin); V_ln <- fit_lin$cov_mat; beta_ln <- beta_ln[rownames(V_ln)]

# covariate column means (shared by all standardization rows)
mk_mm <- function(form, data) model.matrix(form, data)
cov_mean_ns <- colMeans(mk_mm(f_ns, sdat))
cov_mean_ln <- colMeans(mk_mm(f_lin, sdat))
bterms_ns <- grep("^nsEVI|^nsO3", names(beta_ns), value = TRUE)
eterms <- grep("^nsEVI", names(beta_ns), value = TRUE)
oterms <- grep("^nsO3",  names(beta_ns), value = TRUE)

# standardization for one exposure at grid value x (other exposure marginal)
curve_at <- function(xvals, basis, bnames, focal_is_evi) {
  out <- lapply(xvals, function(x) {
    B_new <- predict(basis, rep(x, nrow(sdat)))
    a <- cov_mean_ns
    if (focal_is_evi) a[bnames] <- colMeans(B_new) else a[bnames] <- colMeans(B_new)
    a[bnames] <- colMeans(B_new)
    est <- sum(a * beta_ns)
    se <- sqrt(as.numeric(t(a) %*% V_ns %*% a))
    c(est = est, se = se)
  })
  do.call(rbind, out)
}
lin_at <- function(xvals, term, other_mean) {
  a <- cov_mean_ln
  sapply(xvals, function(x) {
    a[term] <- x
    est <- sum(a * beta_ln)
    se <- sqrt(as.numeric(t(a) %*% V_ln %*% a))
    c(est = est, se = se)
  }) |> t()
}

q_e <- quantile(md$EVIx, c(.01, .25, .5, .75, .99))
q_o <- quantile(md$O3x,  c(.01, .25, .5, .75, .99))
grid_e <- seq(q_e[1], q_e[5], length.out = 120)
grid_o <- seq(q_o[1], q_o[5], length.out = 120)

ce <- curve_at(grid_e, B_e, eterms, TRUE)
co <- curve_at(grid_o, B_o, oterms, FALSE)
le <- lin_at(grid_e, "EVIx"); lo <- lin_at(grid_o, "O3x")

curves <- bind_rows(
  data.frame(exposure = "EVI", x = grid_e, fit = ce[, 1], se = ce[, 2],
             lin = le[, 1]),
  data.frame(exposure = "O3",  x = grid_o, fit = co[, 1], se = co[, 2],
             lin = lo[, 1])) %>%
  mutate(lci = fit - 1.96 * se, uci = fit + 1.96 * se)
write.csv(curves, "output/tables/figure4_exposure_response_v5.csv",
          row.names = FALSE)

# ---- joint contour (linear dual model, marginal standardization) -------------
g1 <- seq(q_e[1], q_e[5], length.out = 70)
g2 <- seq(q_o[1], q_o[5], length.out = 70)
cg <- expand.grid(EVI = g1, O3 = g2)
base_a <- cov_mean_ln
pred_grid <- apply(cg, 1, function(r) {
  a <- base_a; a["EVIx"] <- r[["EVI"]]; a["O3x"] <- r[["O3"]]
  sum(a * beta_ln)
})
cg$fit <- as.numeric(pred_grid)

# mask cells far from the data cloud (standardised nearest-neighbour distance)
z_dat <- scale(md[, c("EVIx", "O3x")])
z_grid <- scale(cg[, c("EVI", "O3")], center = attr(z_dat, "scaled:center"),
                scale = attr(z_dat, "scaled:scale"))
set.seed(42)
nn_idx <- sample.int(nrow(z_dat), min(2000, nrow(z_dat)))
nn_dist <- apply(z_grid, 1, function(g) {
  min(sqrt(colSums((t(z_dat[nn_idx, ]) - g)^2)))
})
cut <- quantile(nn_dist, 0.90)
keep <- nn_dist <= cut
cg$fit[!keep] <- NA
write.csv(cg, "output/tables/figure4_joint_contour_v5.csv", row.names = FALSE)
cat("Contour grid kept:", sum(!is.na(cg$fit)), "of", nrow(cg), "cells\n")

# ---- panels ------------------------------------------------------------------
pal <- list(evi = "#009E73", o3 = "#D55E00")
mk_curve_panel <- function(dsub, col, xlab, tag, q) {
  rng <- range(dsub$uci, dsub$lci, na.rm = TRUE)
  ggplot(dsub, aes(x = x)) +
    annotate("rect", xmin = q[2], xmax = q[4], ymin = -Inf, ymax = Inf,
             fill = "#F0F0F0", alpha = 0.6) +
    geom_rug(data = md, aes(x = if (tag == "a") EVIx else O3x), inherit.aes = FALSE,
             sides = "b", length = unit(0.02, "npc"), alpha = 0.25,
             colour = "#555555", linewidth = 0.2) +
    geom_ribbon(aes(ymin = lci, ymax = uci), fill = col, alpha = 0.18, colour = NA) +
    geom_line(aes(y = lin), linetype = "dashed", colour = "#555555", linewidth = 0.4) +
    geom_line(aes(y = fit), colour = col, linewidth = 0.7) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "#BBBBBB", linewidth = 0.3) +
    labs(x = xlab, y = "Average predicted BAA (years)", tag = tag) +
    coord_cartesian(ylim = rng) +
    theme_nature()
}
p_a <- mk_curve_panel(curves %>% filter(exposure == "EVI"), pal$evi,
                      expression("4-year mean EVI (1 km), IQR band shaded"), "a", q_e)
p_b <- mk_curve_panel(curves %>% filter(exposure == "O3"), pal$o3,
                      expression("4-year mean O"[3]*" (1 km), IQR band shaded"), "b", q_o)

set.seed(7)
pts <- md[sample.int(nrow(md), 1800), ]
p_c <- ggplot() +
  geom_raster(data = cg, aes(x = EVI, y = O3, fill = fit), interpolate = TRUE) +
  geom_contour(data = cg, aes(x = EVI, y = O3, z = fit), colour = "white",
               linewidth = 0.25, alpha = 0.8) +
  geom_point(data = pts, aes(x = EVIx, y = O3x), size = 0.25, alpha = 0.25,
             colour = "black") +
  scale_fill_gradient2(low = "#0072B2", mid = "#F7F7F7", high = "#D55E00",
                       midpoint = 0, name = "Pred. BAA (years)", na.value = "grey92") +
  labs(x = "4-year mean EVI (1 km)",
       y = expression("4-year mean O"[3]*" (1 km, "*mu*"g m"^-3*")"),
       tag = "c") +
  theme_nature() +
  theme(legend.position = "right",
        legend.key.height = unit(4, "mm"), legend.key.width = unit(2.2, "mm"),
        legend.title = element_text(size = 6), legend.text = element_text(size = 5.5))

fig4 <- (p_a | p_b | p_c) + plot_layout(widths = c(1, 1, 1.15))

save_all <- function(plot, stem, width_mm = 183, height_mm = 78) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h,
                       family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = 600, compression = "lzw"); print(plot); dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = w, height = h, units = "in",
                res = 300); print(plot); dev.off()
}
save_all(fig4, "output/figures/Figure4_exposure_response")
cat("\nDone. output/figures/Figure4_exposure_response.{svg,pdf,tiff,png}\n")

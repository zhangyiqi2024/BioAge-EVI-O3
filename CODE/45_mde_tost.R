# =============================================================================
# 45_mde_tost.R  (v5 revision, task T7)
# Quantifies the null difference-in-differences (joint 4-group interaction):
#   - MDE: minimum detectable interaction at 80% power, alpha=0.05 (2.8 x SE)
#   - TOST: two one-sided tests against an equivalence bound of +/- 1 year
# Input: output/tables/joint_association_contrasts_v4.csv (v4 authoritative)
# =============================================================================

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/45_mde_tost.log", split = TRUE)
on.exit(sink(), add = TRUE)

ct <- read.csv("output/tables/joint_association_contrasts_v4.csv",
               stringsAsFactors = FALSE)
did <- ct[ct$contrast == "Difference-in-differences interaction", ]
stopifnot(nrow(did) == 2)

BOUND <- 1.0  # equivalence bound in years of BAA
rows <- list()
for (i in seq_len(nrow(did))) {
  est <- did$estimate[i]; se <- did$se[i]
  mde <- 2.8 * se
  # TOST: H0a: est <= -BOUND  vs H1a: est > -BOUND
  p_lower <- pnorm((est + BOUND) / se, lower.tail = FALSE)  # reject if est far above -bound
  # H0b: est >= +BOUND vs H1b: est < +BOUND
  p_upper <- pnorm((est - BOUND) / se, lower.tail = TRUE)
  p_tost <- max(p_lower, p_upper)
  equivalent <- p_tost < 0.05
  rows[[i]] <- data.frame(
    scale = did$scale[i], estimate = est, se = se,
    lci = est - 1.96 * se, uci = est + 1.96 * se,
    mde_80pct = mde, bound = BOUND,
    tost_p = p_tost, equivalent_at_bound = equivalent,
    interpretation = sprintf(
      "Data rule out interactions larger than %.2f years (MDE, 80%% power); %s",
      mde, ifelse(equivalent,
                  sprintf("equivalence within +/- %.1f year is established (TOST P=%.4f)", BOUND, p_tost),
                  sprintf("equivalence within +/- %.1f year cannot be established (TOST P=%.4f)", BOUND, p_tost))))
  cat(rows[[i]]$scale, ": DiD =", round(est, 3), "SE =", round(se, 3),
      "| MDE =", round(mde, 3), "| TOST p =", round(p_tost, 4),
      "| equivalent:", equivalent, "\n")
}
out <- do.call(rbind, rows)
write.csv(out, "output/tables/did_mde_tost_v5.csv", row.names = FALSE)
cat("\nDone. output/tables/did_mde_tost_v5.csv\n")

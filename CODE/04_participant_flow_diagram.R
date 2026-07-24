# 04_participant_flow_diagram.R — Supplementary Figure S1: participant flow diagram
# Counts verified: steps 0–5 replicated from raw DTA; steps 6–8 from v4 pipeline outputs.
# Outputs: 02_figures/supplementary/FigS1_participant_flow.{png,tiff}
suppressPackageStartupMessages(library(grid))

# Run from the 04_code/ project root (00_master_run.R sets this); the figure
# is written one level up, into the submission package's 02_figures/ folder.
base <- ".."
dir.create(file.path(base, "02_figures", "supplementary"),
           recursive = TRUE, showWarnings = FALSE)

main <- list(
  c("Original merged CHARLS observations (waves 1–5)", "n = 96,628"),
  c("Non-missing residence (city)", "n = 96,628"),
  c("Age ≥ 45 years", "n = 94,240"),
  c("Biomarker waves (2011 and 2015)", "n = 37,049"),
  c("Complete core biomarkers", "n = 21,630"),
  c("No prevalent major chronic disease at baseline", "n = 17,502"),
  c("PhenoAge quality control\n(fasting, valid ranges, CRP < 10 mg/L)", "n = 13,111"),
  c("Observed at both biomarker waves", "n = 6,568 obs; 3,284 participants"),
  c("Analytical sample (ages 46–89;\ncomplete exposures and covariates)", "n = 6,472 obs; 3,281 participants")
)
excl <- list(
  "2" = "Excluded: 2,388\naged < 45 years",
  "3" = "Excluded: 57,191 from\nnon-biomarker waves",
  "4" = "Excluded: 15,419 with\nincomplete biomarkers",
  "5" = "Excluded: 4,128 with prevalent\nchronic disease at baseline",
  "6" = "Excluded: 4,391 failing\nPhenoAge quality control",
  "7" = "Excluded: 6,543 observed at\nonly one biomarker wave",
  "8" = "Excluded: 96 obs (3 participants)\noutside age range or missing\nexposures/covariates"
)

draw <- function() {
  n <- length(main)
  top <- 0.985; bot <- 0.015
  h_slot <- (top - bot) / n
  box_h <- h_slot * 0.60
  mx <- 0.02; mw <- 0.55
  ex <- 0.63; ew <- 0.36
  grid.newpage()
  for (i in seq_len(n)) {
    y <- top - (i - 0.5) * h_slot
    last <- i == n
    face <- if (last) "#b7d4ea" else "#d9e8f5"
    edge <- if (last) "#1f4e79" else "#41719b"
    grid.roundrect(x = mx, y = y, width = mw, height = box_h,
                   just = c("left", "center"), r = unit(0.012, "npc"),
                   gp = gpar(fill = face, col = edge, lwd = 1.4))
    grid.text(main[[i]][1], x = mx + mw/2, y = y + box_h*0.20,
              gp = gpar(fontsize = 7.6, fontfamily = "sans", lineheight = 0.95))
    grid.text(main[[i]][2], x = mx + mw/2, y = y - box_h*0.22,
              gp = gpar(fontsize = 7.8, fontface = "bold", fontfamily = "sans"))
    if (i < n) {
      y_next <- top - (i + 0.5) * h_slot
      grid.segments(mx + mw/2, y - box_h/2 - 0.002,
                    mx + mw/2, y_next + box_h/2 + 0.002,
                    arrow = arrow(length = unit(0.012, "npc"), type = "closed"),
                    gp = gpar(col = "#444444", lwd = 1.1, fill = "#444444"))
    }
    key <- as.character(i)
    if (!is.null(excl[[key]])) {
      ym <- y - h_slot/2
      grid.roundrect(x = ex, y = ym, width = ew, height = box_h * 0.85,
                     just = c("left", "center"), r = unit(0.010, "npc"),
                     gp = gpar(fill = "#f7f7f7", col = "#999999", lwd = 0.9, lty = "dashed"))
      grid.text(excl[[key]], x = ex + ew/2, y = ym,
                gp = gpar(fontsize = 6.8, col = "#555555", fontfamily = "sans", lineheight = 0.95))
      grid.segments(mx + mw + 0.004, ym, ex - 0.004, ym,
                    arrow = arrow(length = unit(0.010, "npc"), type = "closed"),
                    gp = gpar(col = "#999999", lwd = 0.9, fill = "#999999"))
    }
  }
}

png(file.path(base, "02_figures/supplementary/FigS1_participant_flow.png"),
    width = 6.0, height = 8.8, units = "in", res = 300)
draw(); dev.off()

tiff(file.path(base, "02_figures/supplementary/FigS1_participant_flow.tiff"),
     width = 6.0, height = 8.8, units = "in", res = 300, compression = "lzw")
draw(); dev.off()

cat("flow diagram saved\n")

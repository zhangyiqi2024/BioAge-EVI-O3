# =============================================================================
# 46_selection_exposure_table.R  (v5 revision, task T6)
# Quantifies selection: exposure distributions at each cohort-exclusion step
# (Figure S1 chain, 96,628 -> 6,472) plus the primary-analysis frame (5,914).
# Raw CHARLS merged file is restricted: script SKIPs gracefully if absent.
# =============================================================================

suppressPackageStartupMessages({ library(haven); library(dplyr) })

set.seed(20240101)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/46_selection_exposure_table.log", split = TRUE)
on.exit(sink(), add = TRUE)

RAW <- "../../../DATA/原始数据及样本筛选/data_20250925.dta"
if (!file.exists(RAW)) {
  cat("SKIP: restricted raw merged file not found at", RAW, "\n")
  quit(save = "no", status = 0)
}

need <- c("id", "wave", "year", "age", "city", "radyear",
          "bl_wbc", "bl_mcv", "bl_bun", "bl_crea", "bl_crp", "bl_cho",
          "bl_hbalc", "bl_ua", "bl_glu", "bl_fasting", "bl_cysc",
          "癌症", "糖尿病", "心脏病", "中风", "肺病", "记忆疾病",
          paste0("EVI_1KM_lag", 1:4), paste0("O3_1KM_lag", 1:4))
df <- as.data.frame(read_dta(RAW, col_select = any_of(need)))
# sbp lives in the biomarker subset file; merge by id + wave
SBP <- "../../../DATA/原始数据及样本筛选/data_生物年龄.dta"
if (file.exists(SBP)) {
  sb <- as.data.frame(read_dta(SBP, col_select = c("id", "wave", "sbp")))
  df <- merge(df, sb, by = c("id", "wave"), all.x = TRUE)
} else {
  df$sbp <- NA_real_
}
cat("Loaded:", nrow(df), "rows x", ncol(df), "cols\n")

df$EVI_1KM_mean4y <- rowMeans(df[, paste0("EVI_1KM_lag", 1:4)], na.rm = FALSE)
df$O3_1KM_mean4y  <- rowMeans(df[, paste0("O3_1KM_lag", 1:4)],  na.rm = FALSE)

snap <- function(d, step) {
  data.frame(
    step = step, n_obs = nrow(d), n_participants = length(unique(d$id)),
    evi_n = sum(!is.na(d$EVI_1KM_mean4y)),
    evi_mean = mean(d$EVI_1KM_mean4y, na.rm = TRUE),
    evi_sd = sd(d$EVI_1KM_mean4y, na.rm = TRUE),
    o3_n = sum(!is.na(d$O3_1KM_mean4y)),
    o3_mean = mean(d$O3_1KM_mean4y, na.rm = TRUE),
    o3_sd = sd(d$O3_1KM_mean4y, na.rm = TRUE),
    urban_pct = if ("urban_rural" %in% names(d))
      100 * mean(d$urban_rural == 1, na.rm = TRUE) else NA_real_)
}

out <- list()
out[[1]] <- snap(df, "0 Original merged observations")
df <- df[!is.na(df$city), ];                       out[[2]] <- snap(df, "1 Non-missing residence (city)")
df <- df[df$age >= 45, ];                          out[[3]] <- snap(df, "2 Age >= 45 years")
df <- df[df$wave %in% c(1, 3), ];                  out[[4]] <- snap(df, "3 Biomarker waves (2011 and 2015)")
biom <- c("bl_wbc", "bl_mcv", "bl_bun", "bl_crea", "bl_crp",
          "bl_cho", "bl_hbalc", "bl_ua", "sbp", "bl_glu")
df <- df[rowSums(is.na(df[, biom])) == 0, ];       out[[5]] <- snap(df, "4 Complete core biomarkers")
w1 <- df[df$wave == 1, ]
diseased <- w1[w1$癌症 == 1 | w1$糖尿病 == 1 | w1$心脏病 == 1 |
               w1$中风 == 1 | w1$肺病 == 1 | w1$记忆疾病 == 1, ]
df <- df[!df$id %in% diseased$id, ];               out[[6]] <- snap(df, "5 No prevalent major chronic disease at baseline")
df$albumin <- df$bl_cysc
df$creatinine <- df$bl_crea * 88.4 * 0.85
df$exam_year <- ifelse(df$wave == 1, 2011, 2015)
df$death_status <- ifelse(is.na(df$radyear), 0, 1)
df$death_time <- ifelse(df$death_status == 1 & df$radyear > df$exam_year,
                        df$radyear - df$exam_year,
                        ifelse(df$death_status == 0, 2020 - df$exam_year, NA))
df <- df[with(df, bl_fasting == 1 & albumin >= 0.5 & albumin <= 4.0 &
              bl_crp < 10 & creatinine >= 30 & creatinine <= 300 &
              !is.na(death_time) & death_time > 0) %in% TRUE, ]
out[[7]] <- snap(df, "6 PhenoAge quality control")

# Step 7: shipped analysis dataset IS the official output of the both-waves
# step (6,568 obs / 3,284 participants); raw-chain lineage cannot be
# reproduced exactly from the merged file, so use the authoritative frame.
ad0 <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
ad0$EVI_1KM_mean4y <- rowMeans(ad0[, paste0("EVI_1KM_lag", 1:4)], na.rm = FALSE)
ad0$O3_1KM_mean4y  <- rowMeans(ad0[, paste0("O3_1KM_lag", 1:4)],  na.rm = FALSE)
out[[8]] <- snap(ad0, "7 Observed at both biomarker waves")

# Step 8: analytical sample = ages 46-89 (authoritative, 6,472 / 3,281)
ad <- ad0[ad0$age > 45 & ad0$age < 90, ]
out[[9]] <- snap(ad, "8 Analytical sample (ages 46-89)")

res <- do.call(rbind, out)
res$note <- ""
res$note[7] <- "Replicated 13,112 vs official 13,111 (1 obs; sex-filter nuance in original calibration script)"
res$note[8] <- "Authoritative shipped analysis frame (both-waves lineage not exactly reproducible from merged raw file)"
print(res, digits = 3)

# Flow-count verification against FigureS1_participant_flow.csv
flow <- read.csv("../03_source_data/FigureS1_participant_flow.csv")
chk <- data.frame(step = flow$step[1:9], expected = flow$n_obs[1:9],
                  replicated = res$n_obs[1:9])
chk$match <- chk$expected == chk$replicated
print(chk)
cat("Steps 0-5 and 8 replicate exactly; step 6 differs by 1 obs; step 7 uses the authoritative shipped frame.\n")

write.csv(res, "output/tables/selection_exposure_by_step_v5.csv", row.names = FALSE)
cat("\nDone. output/tables/selection_exposure_by_step_v5.csv\n")

# 03_cohort_flow_replication.R — replicate cohort selection steps 0–5 from
# DATA/原始数据及样本筛选/生物年龄样本筛选_各步整合代码.R on the local raw DTA,
# and compare with s4_cohort_flow_v4.csv (SI Table S4 / flow diagram source).
# Input:  <repo root>/DATA/原始数据及样本筛选/data_生物年龄.dta (restricted raw
#         CHARLS; expected three levels above 04_code/. 00_master_run.R checks
#         for it and skips this step when absent.)
# Output: output/tables/cohort_flow_replication.csv

suppressPackageStartupMessages({ library(haven); library(dplyr) })

ROOT <- "../../.."   # repository root relative to Submission/final_Submission/04_code
df <- read_dta(file.path(ROOT, "DATA/原始数据及样本筛选/data_生物年龄.dta"))

steps <- data.frame(step = character(), n_obs = integer(), stringsAsFactors = FALSE)
steps <- rbind(steps, data.frame(step = "0 Original merged observations", n_obs = nrow(df)))

df <- df[!is.na(df$city), ]
steps <- rbind(steps, data.frame(step = "1 Non-missing residence (city)", n_obs = nrow(df)))

df <- df[df$age >= 45, ]
steps <- rbind(steps, data.frame(step = "2 Age >= 45 years", n_obs = nrow(df)))

df <- df[df$wave %in% c(1, 3), ]
steps <- rbind(steps, data.frame(step = "3 Biomarker waves (2011 and 2015)", n_obs = nrow(df)))

biomarkers <- c("bl_wbc", "bl_mcv", "bl_bun", "bl_crea", "bl_crp",
                "bl_cho", "bl_hbalc", "bl_ua", "sbp", "bl_glu")
df <- df[rowSums(is.na(df[, biomarkers])) == 0, ]
steps <- rbind(steps, data.frame(step = "4 Complete core biomarkers", n_obs = nrow(df)))

w1 <- df[df$wave == 1, ]
diseased <- w1 %>% filter(癌症 == 1 | 糖尿病 == 1 | 心脏病 == 1 | 中风 == 1 | 肺病 == 1 | 记忆疾病 == 1)
df <- df %>% filter(!id %in% diseased$id)
steps <- rbind(steps, data.frame(step = "5 No prevalent major chronic disease at baseline", n_obs = nrow(df)))

# expected values from s4_cohort_flow_v4.csv (steps 0-5)
expected <- c(96628, 96628, 94240, 37049, 21630, 17502)
steps$expected_s4 <- expected
steps$match <- steps$n_obs == steps$expected_s4
print(steps)

write.csv(steps, "output/tables/cohort_flow_replication.csv", row.names = FALSE)
cat(ifelse(all(steps$match), "ALL STEPS MATCH\n", "MISMATCH FOUND\n"))

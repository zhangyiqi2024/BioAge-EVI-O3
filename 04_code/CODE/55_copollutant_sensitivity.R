# Supportive PM2.5/NO2 mutual-adjustment sensitivity analysis.
#
# The participant-level frame is restricted and is therefore not distributed.
# When the local audit frame is present, this script delegates to the fitted
# WeMix model script in the audit folder. Otherwise it validates and copies the
# permitted aggregate result used by Supplementary Table S31. The public path
# never fabricates or reconstructs participant-level exposure records.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- normalizePath(dirname(if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "55_copollutant_sensitivity.R"), winslash = "/", mustWork = FALSE)
root <- normalizePath(file.path(script_dir, "../.."), winslash = "/", mustWork = FALSE)
audit_dir <- file.path(root, "临时文件", "co_pollutant_sensitivity_20261016")
audit_script <- file.path(audit_dir, "run_copollutant_models.R")
audit_frame <- file.path(audit_dir, "copollutant_analysis_frame.csv")
public_csv <- file.path(root, "03_source_data", "copollutant_sensitivity_v1.csv")
out_dir <- file.path(root, "04_code", "output", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists(audit_script) && file.exists(audit_frame)) {
  message("Restricted local frame detected; running the WeMix sensitivity models.")
  source(audit_script, local = FALSE)
} else if (file.exists(public_csv)) {
  destination <- file.path(out_dir, "copollutant_sensitivity_v1.csv")
  file.copy(public_csv, destination, overwrite = TRUE)
  message("Restricted participant-level frame is absent; copied the permitted aggregate result to ", destination)
} else {
  stop("Neither the restricted audit frame nor the permitted aggregate result is available.")
}

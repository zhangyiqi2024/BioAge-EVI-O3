# Install the declared direct dependencies into an existing R environment.
# Usage: Rscript requirements.R
# Exact tested versions are recorded in requirements_r_versions.csv.
manifest <- read.csv("requirements_r_versions.csv", stringsAsFactors = FALSE)
missing <- manifest$package[!vapply(manifest$package, requireNamespace,
                                  logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")

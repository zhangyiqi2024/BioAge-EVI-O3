# Dependency and numerical-environment report; does not install packages.
manifest <- read.csv("requirements_r_versions.csv", stringsAsFactors = FALSE)
installed <- installed.packages()
present <- manifest$package %in% rownames(installed)
if (!all(present)) stop("Missing R packages: ", paste(manifest$package[!present], collapse = ", "))
dir.create("logs", recursive = TRUE, showWarnings = FALSE)
versions <- data.frame(package = manifest$package,
    expected = manifest$version,
    installed = installed[manifest$package, "Version"])
write.csv(versions, "logs/dependency_versions.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "logs/session_info.txt")
different <- versions$expected != versions$installed
if (any(different)) warning("Version differences: ", paste(versions$package[different], collapse = ", "))
cat("R dependency preflight passed for", nrow(versions), "packages.\n")

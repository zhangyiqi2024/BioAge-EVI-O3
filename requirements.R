# R package requirements for reproducible analysis
# Install with: install.packages(pkgs) and BiocManager::install() where noted

pkgs <- c(
  # Data handling
  "haven",       # Read/write Stata .dta files
  "dplyr",       # Data manipulation
  "tidyr",       # Data reshaping
  "readr",       # Read CSV
  "readxl",      # Read Excel
  "writexl",     # Write Excel

  # Visualization
  "ggplot2",     # Core plotting
  "ggpubr",      # Publication-ready plots
  "patchwork",   # Combine multiple ggplots
  "scales",      # Axis formatting
  "RColorBrewer",# Color palettes
  "viridis",     # Colorblind-friendly palettes

  # Statistical modeling
  "lme4",        # Linear mixed-effects models
  "lmerTest",    # P-values for lme4
  "metafor",     # Meta-analysis / forest plots
  "qgcomp",      # Quantile g-computation
  "mgcv",        # GAM / nonlinear modeling
  "splines",     # Spline basis functions

  # Machine learning
  "lightgbm",    # Gradient boosting
  "shapviz",     # SHAP visualization
  "shapr",       # SHAP values

  # Tables & reporting
  "gtsummary",   # Summary tables
  "flextable",   # Formatted tables
  "officer",     # Office document generation
  "kableExtra",  # LaTeX/HTML tables

  # Utilities
  "stringr",     # String manipulation
  "lubridate",   # Date handling
  "broom",       # Tidy model outputs
  "broom.mixed", # Tidy lme4 outputs
  "purrr"        # Functional programming
)

# Install all
install.packages(pkgs)

# Note: lightgbm may need additional system dependencies on Linux
# See: https://lightgbm.readthedocs.io/en/latest/Installation-Guide.html

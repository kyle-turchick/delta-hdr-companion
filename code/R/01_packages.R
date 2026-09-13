# Checks for the required R packages and loads them for the workflows.

# Package setup ####
# ----------------------------- #
required_packages <- c(
  "readr", "dplyr", "tidyr", "tibble", "purrr", "terra", "broom",
  "ggplot2", "patchwork", "hexbin", "scales", "rlang",
  "officer", "flextable"
)
missing_packages <- required_packages[
  !vapply(
    required_packages,
    function(pkg) suppressWarnings(requireNamespace(pkg, quietly = TRUE)),
    logical(1)
  )
]
if (length(missing_packages) > 0) {
  if (!INSTALL_MISSING) {
    stop("Required R packages are missing: ", paste(missing_packages, collapse = ", "))
  }
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
suppressWarnings(suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(terra)
  library(broom)
  library(ggplot2)
  library(patchwork)
  library(hexbin)
  library(scales)
}))

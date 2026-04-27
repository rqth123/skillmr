#!/usr/bin/env Rscript
#
# Installation script for MRAgent R dependencies
# Usage: Rscript install_dependencies.R
#

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/")
  }
}

# Required packages for standard MR
required_packages <- c(
  "TwoSampleMR",
  "ieugwasr",
  "dplyr",
  "jsonlite",
  "pdftools"
)

# Required packages for MRlap (optional)
mrlap_packages <- c(
  "MRlap",
  "httr",
  "vcfR"
)

# Install main dependencies
cat("Installing main dependencies...\n")
for (pkg in required_packages) {
  install_if_missing(pkg)
}

# Install MRlap dependencies if needed
cat("\nDo you want to install MRlap dependencies for sample overlap correction? (y/n): ")
answer <- readLines(n = 1)
if (tolower(answer) == "y") {
  cat("Installing MRlap dependencies...\n")
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")
  }
  BiocManager::install("vcfR")
  for (pkg in mrlap_packages) {
    if (pkg != "vcfR") {
      install_if_missing(pkg)
    }
  }
}

cat("\nInstallation complete!\n")
cat("Required packages installed:\n")
for (pkg in required_packages) {
  status <- ifelse(requireNamespace(pkg, quietly = TRUE), "✓", "✗")
  cat(sprintf("  %s %s\n", status, pkg))
}

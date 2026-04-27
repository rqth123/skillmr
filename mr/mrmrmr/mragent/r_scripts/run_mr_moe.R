#!/usr/bin/env Rscript
#
# Mixture of Experts MR analysis - standalone R script
# Input: Exposure_id, Outcome_id, output_path, gwas_token
# Output: CSV results and PDF plots
#

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_mr_moe.R <Exposure_id> <Outcome_id> <output_path> [gwas_token]\n")
}

Exposure_id <- args[1]
Outcome_id <- args[2]
path <- args[3]
gwas_token <- if (length(args) >= 4) args[4] else Sys.getenv("OPENGWAS_JWT", "")

# Set environment
Sys.setenv(OPENGWAS_JWT = gwas_token)

# Load libraries
if (!requireNamespace("TwoSampleMR", quietly = TRUE) || !requireNamespace("ieugwasr", quietly = TRUE)) {
  stop("Required packages not installed. Please run install_dependencies.R first")
}
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)

# Extract instruments
p_value <- 5e-08
exposure_dat <- extract_instruments(
  outcomes = Exposure_id,
  p1 = p_value,
  clump = TRUE,
  r2 = 0.001,
  kb = 5000
)

num_rows <- nrow(exposure_dat)
message(sprintf("%d SNPs extracted", num_rows))

# Relax p-value if not enough SNPs
if (num_rows < 5) {
  p_value <- 5e-06
  exposure_dat <- extract_instruments(
    outcomes = Exposure_id,
    p1 = p_value,
    clump = TRUE,
    r2 = 0.001,
    kb = 5000
  )
  num_rows <- nrow(exposure_dat)
  message(sprintf("Retried with p < 5e-06: %d SNPs extracted", num_rows))
}

if (num_rows < 3) {
  stop(sprintf("Failed to extract instruments. Only %d SNPs found. Need at least 3.", num_rows))
}

# Extract outcome data
message("Extracting outcome data...")
outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = Outcome_id)
dat <- harmonise_data(exposure_dat, outcome_dat)
outTab <- dat[dat$mr_keep == TRUE, ]

# Validate output
if (nrow(outTab) < 3) {
  stop(sprintf("After harmonisation, only %d SNPs remain. Need at least 3.", nrow(outTab)))
}

# Write SNP data
snp_file <- file.path(path, "table.SNP.csv")
write.csv(outTab, file = snp_file, row.names = FALSE)
message(sprintf("Saved SNP data to %s", snp_file))

# ====== MR Mixture of Experts analysis ======
message("Running MR Mixture of Experts analysis...")

# Apply all MR methods
r <- mr_wrapper(dat)

# Load the rf object containing the trained models
# Note: rf.rdata should be in the working directory
if (!file.exists("rf.rdata")) {
  stop("rf.rdata not found. This is required for MR-MOE analysis.")
}
load("rf.rdata")

# Update the results with mixture of experts
r <- mr_moe(r, rf)

# Save results
mr_file <- file.path(path, "MR.MRresult.csv")
write.csv(r[[1]]$estimates, file = mr_file, row.names = FALSE)

heter_file <- file.path(path, "MR.heterogeneity.csv")
write.csv(r[[1]]$heterogeneity, file = heter_file, row.names = FALSE)

pleio_file <- file.path(path, "MR.table.pleiotropy.csv")
write.csv(pleioTab, file = pleio_file, row.names = FALSE)

# Heterogeneity analysis
heterTab <- mr_heterogeneity(dat)

# Pleiotropy test
pleioTab <- mr_pleiotropy_test(dat)

# ====== Generate plots ======
message("Generating plots...")

# Scatter plot
pdf(file = file.path(path, "pic.scatter_plot.pdf"), width = 7.5, height = 7)
mr_scatter_plot(r[[1]]$estimates, dat)
dev.off()

# Forest plot
res_single <- mr_singlesnp(dat)
pdf(file = file.path(path, "pic.forest.pdf"), width = 7, height = 5.5)
mr_forest_plot(res_single)
dev.off()

# Funnel plot
pdf(file = file.path(path, "pic.funnel_plot.pdf"), width = 7, height = 6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

# Leave-one-out sensitivity analysis
pdf(file = file.path(path, "pic.leaveoneout.pdf"), width = 7, height = 5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()

message("MR-MOE analysis completed successfully!")
message(sprintf("Output directory: %s", path))
message(sprintf("Final SNP count: %d", nrow(outTab)))

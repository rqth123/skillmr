#!/usr/bin/env Rscript
#
# MRlap analysis for sample overlap correction - standalone R script
# Input: Exposure_id, Outcome_id, output_path, N_exposure, N_outcome
# Output: JSON results
#

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript run_mrlap.R <Exposure_id> <Outcome_id> <output_path> <N_exposure> <N_outcome>\n")
}

Exposure_id <- args[1]
Outcome_id <- args[2]
path <- args[3]
N_exposure <- as.numeric(args[4])
N_outcome <- as.numeric(args[5])

# Check required packages
required_pkgs <- c("MRlap", "httr", "vcfR", "jsonlite")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' not installed. Please install it first.", pkg))
  }
}

library(httr)
library(vcfR)
library(MRlap)
library(jsonlite)

# Download function
download_vcf <- function(base_url, file_name, dest_dir = ".", min_size = 1 * 1024 * 1024) {
  file_id <- sub("(.*)\\.vcf\\.gz", "\\1", file_name)
  file_url <- paste0(base_url, file_id, "/", file_name)
  message(sprintf("Attempting to download from: %s", file_url))

  dest_file <- file.path(dest_dir, file_name)
  message(sprintf("Saving to: %s", dest_file))

  if (file.exists(dest_file)) {
    file_info <- file.info(dest_file)
    file_size <- file_info$size
    if (file_size >= min_size) {
      message(sprintf("File already exists and is larger than %.2f MB. Skipping download.",
                      round(min_size / (1024 * 1024), 2)))
      return(TRUE)
    } else {
      message(sprintf("File exists but is too small (%.2f MB). Redownloading...",
                      round(file_size / (1024 * 1024), 2)))
    }
  }

  tryCatch({
    download.file(file_url, destfile = dest_file, method = "curl", mode = "wb")
    file_info <- file.info(dest_file)
    file_size <- file_info$size
    if (file_size < min_size) {
      message(sprintf("File size is too small (%.2f MB). Deleting file: %s",
                      round(file_size / (1024 * 1024), 2), dest_file))
      file.remove(dest_file)
      return(FALSE)
    }
    message(sprintf("Download complete: %s, File size: %.2f MB",
                    dest_file, round(file_size / (1024 * 1024), 2)))
    return(TRUE)
  }, error = function(e) {
    if (grepl("404", e$message)) {
      message(sprintf("File not found (404): %s", file_name))
    } else {
      message(sprintf("Error downloading file: %s, Error message: %s", file_name, e$message))
    }
    return(FALSE)
  })
}

# Extract data from VCF
extract_data_from_vcf <- function(vcf_file, N) {
  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  sample_columns <- colnames(vcf@gt)
  sample_name <- sample_columns[length(sample_columns)]
  message("Available sample columns in the VCF:")
  print(sample_name)

  fix_data <- as.data.frame(vcf@fix)
  fix_data$POS <- as.numeric(fix_data$POS)

  gt_data <- as.data.frame(vcf@gt)
  gt_parsed <- do.call(rbind, strsplit(gt_data[, sample_name], ":"))

  if (nrow(gt_parsed) != nrow(fix_data)) {
    stop("The number of rows in the genotype data does not match the fixed fields. Check the VCF file.")
  }

  df <- data.frame(
    chr = fix_data$CHROM,
    pos = fix_data$POS,
    rsid = fix_data$ID,
    ref = fix_data$REF,
    alt = fix_data$ALT,
    beta = as.numeric(gt_parsed[, 1]),
    se = as.numeric(gt_parsed[, 2]),
    N = N
  )
  df$zscore <- df$beta / df$se

  message("Data extraction complete.")
  print(head(df))
  return(df)
}

# Main analysis
base_url <- "https://gwas.mrcieu.ac.uk/files/"
exposure_file_name <- paste0(Exposure_id, ".vcf.gz")
outcome_file_name <- paste0(Outcome_id, ".vcf.gz")

# LD reference files (should be in working directory or adjust path)
ld_file <- Sys.getenv("MRAGENT_LD_PATH", "./eur_w_ld_chr")
hm3_file <- Sys.getenv("MRAGENT_HM3_PATH", "./w_hm3.snplist")

message(sprintf("Using LD files: %s, %s", ld_file, hm3_file))

# Download and process exposure data
if (!download_vcf(base_url, exposure_file_name)) {
  stop("Failed to download exposure data.")
}
exposure_data <- extract_data_from_vcf(exposure_file_name, N_exposure)

# Download and process outcome data
if (!download_vcf(base_url, outcome_file_name)) {
  stop("Failed to download outcome data.")
}
outcome_data <- extract_data_from_vcf(outcome_file_name, N_outcome)

# Run MRlap
message("Running MRlap analysis...")
result <- MRlap(
  exposure = exposure_data,
  exposure_name = "Exposure",
  outcome = outcome_data,
  outcome_name = "Outcome",
  ld = ld_file,
  hm3 = hm3_file,
  MR_threshold = 5e-6
)

print(result)

# Save results
result_json <- toJSON(result, pretty = TRUE, auto_unbox = TRUE)
output_file <- file.path(path, "MRlap_results.json")
writeLines(result_json, output_file)

message(sprintf("MRlap analysis completed. Results saved to: %s", output_file))

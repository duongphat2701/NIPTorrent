#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(NIPTeR)
  library(optparse)
})

# ----------------------------
# Chi-square variation reduction (χ2VR) for CONTROL GROUP ONLY
# Implements the method described in the NIPTeR manual:
# - For each bin across autosomes, compute chi-square across controls
# - Convert to ~N(0,1): (chi - df)/sqrt(2*df)
# - If above cutoff, divide reads by (chi/df)
# Keeps ALL controls (no removals).
# ----------------------------
chi_correct_controlgroup <- function(nipt_control_group, chi_cutoff = 3.5) {
  cg <- nipt_control_group
  samples <- cg$Samples
  n <- length(samples)
  if (n < 2) stop("Need at least 2 samples for chi correction")

  # Require non strand-separated data: autosomal_chromosome_reads is a list of length 1
  k <- length(samples[[1]]$autosomal_chromosome_reads)
  if (k != 1) {
    stop("chi_correct_controlgroup requires separate_strands = FALSE (single-track autosomal_chromosome_reads).")
  }

  # bins count from first sample
  nbins <- ncol(samples[[1]]$autosomal_chromosome_reads[[1]])
  df <- n - 1

  # autosomes only
  for (chr in 1:22) {
    # Build matrix (n_samples x nbins)
    M <- matrix(0, nrow = n, ncol = nbins)
    for (i in seq_len(n)) {
      M[i, ] <- samples[[i]]$autosomal_chromosome_reads[[1]][chr, ]
    }

    mu <- colMeans(M)
    valid <- mu > 0

    chi <- rep(0, nbins)
    if (any(valid)) {
      mu_mat <- matrix(mu[valid], n, sum(valid), byrow = TRUE)
      chi[valid] <- colSums((M[, valid, drop = FALSE] - mu_mat)^2 / mu_mat)
    }

    # Normalized chi-square score ~ N(0,1)
    z <- rep(-Inf, nbins)
    if (any(valid)) {
      z[valid] <- (chi[valid] - df) / sqrt(2 * df)
    }

    idx <- which(z > chi_cutoff)
    if (length(idx)) {
      factor <- chi[idx] / df
      # apply correction to all samples
      M[, idx] <- M[, idx, drop = FALSE] / matrix(factor, n, length(idx), byrow = TRUE)

      # write back
      for (i in seq_len(n)) {
        samples[[i]]$autosomal_chromosome_reads[[1]][chr, ] <- M[i, ]
      }
    }
  }

  # Update correction status (optional but helpful)
  for (i in seq_len(n)) {
    samples[[i]]$correction_status_autosomal_chromosomes <-
      paste(unique(c(samples[[i]]$correction_status_autosomal_chromosomes, "Chi Corrected")),
            collapse = " + ")
  }
  cg$Samples <- samples
  cg$Correction_status <- paste(unique(c(cg$Correction_status, "Chi Corrected")), collapse = " + ")
  cg
}

# ----------------------------
# Sample name setter to prevent overwriting in as_control_group()
# ----------------------------
set_sample_name_safe <- function(s, nm) {
  # S4 (typical for NIPTeR)
  if (isS4(s)) {
    if ("sample_name" %in% slotNames(s)) {
      slot(s, "sample_name") <- nm
      return(s)
    }
    if ("samplename" %in% slotNames(s)) {
      slot(s, "samplename") <- nm
      return(s)
    }
  }
  # fallback for list-like
  s$sample_name <- nm
  s
}

# ----------------------------
# CLI
# ----------------------------
option_list <- list(
  make_option(c("-l", "--bam_list"), type = "character",
              help = "Text file with paths to BAMs, one per line"),
  make_option(c("-d", "--bam_dir"), type = "character",
              help = "Directory containing BAMs (alternative to --bam_list)"),
  make_option(c("-o", "--out_rds"), type = "character", default = "controlgroup.gc_loess_chi.rds",
              help = "Output RDS path [default: %default]"),
  make_option(c("--diagnostics_tsv"), type = "character", default = "controlgroup.diagnostics.tsv",
              help = "Write diagnose_control_group() aberrant scores to TSV [default: %default]"),
  make_option(c("--chi_cutoff"), type = "double", default = 3.5,
              help = "Chi cutoff for χ2VR [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$bam_list) && is.null(opt$bam_dir)) {
  stop("You must provide either --bam_list or --bam_dir")
}

# ----------------------------
# 1) Collect BAM paths
# ----------------------------
if (!is.null(opt$bam_list)) {
  bam_filepaths <- scan(opt$bam_list, what = "character", quiet = TRUE)
} else {
  bam_filepaths <- list.files(
    path       = opt$bam_dir,
    pattern    = "\\.bam$",
    full.names = TRUE
  )
}

bam_filepaths <- bam_filepaths[file.exists(bam_filepaths)]
if (!length(bam_filepaths)) stop("No .bam files found.")
cat("Found", length(bam_filepaths), "BAM files\n")

# ----------------------------
# 2) Bin BAMs -> NIPTSamples
#    IMPORTANT: separate_strands=FALSE for stable control-group χ2VR build
# ----------------------------
samples <- lapply(bam_filepaths, function(bam) {
  cat("Loading Bam:", bam, "\n")
  
  # Load BAM and check that it's properly processed
  s <- tryCatch(
    bin_bam_sample(bam, do_sort = FALSE, separate_strands = FALSE),
    error = function(e) {
      cat("Error loading BAM:", bam, "Error message:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  # If no sample returned, skip this BAM
  if (is.null(s)) return(NULL)
  
  # Generate sample name based on filename
  nm <- sub("\\.bam$", "", basename(bam))  # Remove .bam extension
  s <- set_sample_name_safe(s, nm)  # Force setting the sample_name
  
  return(s)
})

# Remove any NULL entries that might have occurred due to errors
samples <- samples[!sapply(samples, is.null)]

cat("Control group contains:", length(samples), "samples\n")
sample_names <- vapply(samples, function(s) s$sample_name, character(1))
cat("Sample names:", paste(sample_names, collapse = ", "), "\n")

# Now create the control group
if (length(samples) > 0) {
  control_group <- as_control_group(nipt_samples = samples)
  cat("Control group contains:", length(control_group$Samples), "samples\n")
  cat("Sample names in control group:", paste(names(control_group$Samples), collapse = ", "), "\n")
} else {
  cat("No valid samples found. Check the BAM files for errors.\n")
}

# ----------------------------
# 3) GC correction (fixed: LOESS, exclude sex chromosomes)
# ----------------------------
gc_control_group <- gc_correct(
  nipt_object = control_group,
  method      = "LOESS",
  include_XY  = FALSE
)

# ----------------------------
# 4) Chi-squared variation reduction (fixed TRUE) - CONTROL GROUP ONLY
# ----------------------------
final_control_group <- chi_correct_controlgroup(
  nipt_control_group = gc_control_group,
  chi_cutoff         = opt$chi_cutoff
)

# ----------------------------
# 5) Diagnose (report only; keep all controls)
# ----------------------------
diag <- diagnose_control_group(nipt_control_group = final_control_group)

aberrant_df <- NULL
if (!is.null(diag$abberant_scores)) aberrant_df <- diag$abberant_scores
if (!is.null(diag$aberrant_scores)) aberrant_df <- diag$aberrant_scores

if (!is.null(aberrant_df)) {
  write.table(
    aberrant_df,
    file      = opt$diagnostics_tsv,
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
  )
  cat("Diagnostics written to:", opt$diagnostics_tsv, "\n")
  cat("NOTE: Keeping all controls (no removals).\n")
} else {
  cat("diagnose_control_group() returned no aberrant score table to export.\n")
}

# ----------------------------
# 6) Save
# ----------------------------
saveRDS(final_control_group, opt$out_rds)
cat("Saved control group to:", opt$out_rds, "\n")

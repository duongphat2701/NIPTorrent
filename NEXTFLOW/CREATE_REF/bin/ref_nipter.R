#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(NIPTeR)
  library(optparse)
})

option_list <- list(
  make_option(c("-l", "--bam_list"), type = "character",
              help = "Text file with paths to control BAMs, one per line"),
  make_option(c("-d", "--bam_dir"), type = "character",
              help = "Directory with control BAMs (alternative to --bam_list)"),
  make_option(c("-o", "--out_rds"), type = "character", default = "controlgroup.rds",
              help = "Output RDS path [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$bam_list) && is.null(opt$bam_dir)) {
  stop("You must provide either --bam_list or --bam_dir")
}

## 1) Collect BAM paths
if (!is.null(opt$bam_list)) {
  bam_filepaths <- scan(opt$bam_list, what = "character", quiet = TRUE)
} else {
  bam_filepaths <- list.files(
    path       = opt$bam_dir,
    pattern    = "\\.bam$",
    full.names = TRUE
  )
}

if (!length(bam_filepaths)) {
  stop("No .bam files found.")
}

cat("Found", length(bam_filepaths), "BAM files\n")

## 2) Load BAMs -> NIPTSample -> NIPTControlGroup
samples <- lapply(bam_filepaths, function(bam) {
  cat("Loading Bam:", bam, "\n")
  
  # Load BAM and check that it's properly processed
  s <- bin_bam_sample(
    bam,
    do_sort          = FALSE,
    separate_strands = FALSE
  )
  
# Set sample name if not assigned from BAM (fallback to filename)
set_sample_name_safe <- function(s, nm) {
  if (is.null(s$sample_name)) {
    cat("Warning: No sample_name assigned for BAM. Using fallback:", nm, "\n")
    s$sample_name <- nm  # Force sample name if missing
  }
  return(s)
}

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


# Check how many samples were created and print their names
cat("Control group contains:", length(samples), "samples\n")
sample_names <- vapply(samples, function(s) s$sample_name, character(1))
cat("Sample names:", paste(sample_names, collapse = ", "), "\n")

# Now create the control group and check its contents
control_group <- as_control_group(nipt_samples = samples)
cat("Control group contains:", length(control_group$Samples), "samples\n")
cat("Sample names in control group:", paste(names(control_group$Samples), collapse = ", "), "\n")


## 3) GC correction
gc_control_group <- gc_correct(
  nipt_object = control_group,
  method      = "LOESS",
  include_XY  = FALSE
)

## 4) Diagnose control group
control_group_diagnostics <- diagnose_control_group(
  nipt_control_group = gc_control_group
)

## 5) Remove aberrant controls
aberrant_sample_names <- unique(
  control_group_diagnostics$abberant_scores$sample_name
)

if (length(aberrant_sample_names)) {
  cat("Removing", length(aberrant_sample_names), "aberrant controls:\n")
  print(aberrant_sample_names)
  for (nm in aberrant_sample_names) {
    gc_control_group <- remove_sample_controlgroup(
      samplename         = nm,
      nipt_control_group = gc_control_group
    )
  }
} else {
  cat("No aberrant controls detected\n")
}

## 6) Save
saveRDS(object = gc_control_group, file = opt$out_rds)
cat("Saved control group to:", opt$out_rds, "\n")

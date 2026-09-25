#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(NIPTeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage:\n",
    "  Rscript nipter_allchr.R <controlgroup_rds> <sample_bam> <out_prefix> [chroms]\n\n",
    "chroms (optional): e.g. '1-22' (default), '13,18,21', '1:10'\n\n",
    "Example:\n",
    "  Rscript nipter_allchr.R NIPTeR_control.rds 24PC1151_021.bam 24PC1151_021 1-22\n"
  )
}

control_rds <- args[1]
sample_bam  <- args[2]
out_prefix  <- args[3]
chrom_arg   <- if (length(args) >= 4) args[4] else "1-22"

if (!file.exists(control_rds)) stop("Control group RDS not found: ", control_rds)
if (!file.exists(sample_bam))  stop("Sample BAM not found: ", sample_bam)

# ---------------------------
# helpers
# ---------------------------

parse_chroms <- function(x) {
  x <- gsub("\\s+", "", x)
  if (grepl("^[0-9]+-[0-9]+$", x)) {
    a <- as.integer(strsplit(x, "-", fixed = TRUE)[[1]])
    return(seq.int(a[1], a[2]))
  }
  if (grepl("^[0-9]+:[0-9]+$", x)) {
    a <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    return(seq.int(a[1], a[2]))
  }
  if (grepl("^[0-9]+(,[0-9]+)+$", x) || grepl("^[0-9]+$", x)) {
    return(as.integer(strsplit(x, ",", fixed = TRUE)[[1]]))
  }
  stop("Could not parse chroms: ", x, "  (use '1-22', '1:22', or '1,2,3')")
}

safe_try <- function(expr) {
  tryCatch(expr, error = function(e) e)
}

has_string <- function(x, pattern) {
  if (is.null(x) || !is.character(x)) return(FALSE)
  any(grepl(pattern, x, ignore.case = TRUE))
}

get_z_score <- function(zobj) {
  if (inherits(zobj, "error") || is.null(zobj)) return(NA_real_)
  if (!is.null(zobj$sample_Zscore) && is.numeric(zobj$sample_Zscore) && length(zobj$sample_Zscore) == 1)
    return(zobj$sample_Zscore)
  NA_real_
}

get_ncv_score <- function(ncvobj) {
  if (inherits(ncvobj, "error") || is.null(ncvobj)) return(NA_real_)
  # NCVResult is basically NCVTemplate + sample score appended (name varies across versions/usages).
  # Try the most likely fields first.
  for (nm in c("sample_Zscore", "sample_NCVscore", "sample_NCV_score", "NCV_score_sample", "sample_score")) {
    if (!is.null(ncvobj[[nm]]) && is.numeric(ncvobj[[nm]]) && length(ncvobj[[nm]]) == 1)
      return(ncvobj[[nm]])
  }
  # fallback: first numeric scalar that is NOT focus_chromosome
  nms <- names(ncvobj)
  for (nm in nms) {
    if (nm %in% c("focus_chromosome")) next
    val <- ncvobj[[nm]]
    if (is.numeric(val) && length(val) == 1) return(val)
  }
  NA_real_
}

get_rbz_sets <- function(regobj) {
  if (inherits(regobj, "error") || is.null(regobj)) return(numeric(0))
  ps <- regobj$prediction_statistics
  if (is.null(ps) || is.null(rownames(ps)) || !("Z_score_sample" %in% rownames(ps))) return(numeric(0))
  zvec <- suppressWarnings(as.numeric(ps["Z_score_sample", , drop = TRUE]))
  zvec[is.finite(zvec)]
}

get_rbz_score <- function(regobj, agg = c("mean", "median", "first")) {
  agg <- match.arg(agg)
  zvec <- get_rbz_sets(regobj)
  if (!length(zvec)) return(NA_real_)
  switch(agg,
         mean   = mean(zvec),
         median = median(zvec),
         first  = zvec[1])
}

# ---------------------------
# 1) Load control group
# ---------------------------

control_group <- readRDS(control_rds)

if (is.null(control_group$samples) || length(control_group$samples) < 1) {
  stop("control_group$samples is empty or missing; is this a valid NIPTControlGroup RDS?")
}

# Detect whether controls were binned with separate_strands
ctrl_sep <- FALSE
first_sample <- control_group$samples[[1]]
if (!is.null(first_sample$autosomal_chromosome_reads)) {
  ctrl_sep <- length(first_sample$autosomal_chromosome_reads) == 2
}

# ---------------------------
# 2) Load / bin sample in the SAME mode as controls
# ---------------------------

sample_of_interest <- bin_bam_sample(
  bam_filepath     = sample_bam,
  do_sort          = FALSE,
  separate_strands = ctrl_sep
)

# ---------------------------
# 3) GC correction + match + chi (guard against double correction)
# ---------------------------

cg_status <- control_group$correction_status_autosomal_chromosomes
sm_status <- sample_of_interest$correction_status_autosomal_chromosomes

if (!has_string(cg_status, "GC")) {
  control_group <- gc_correct(nipt_object = control_group, method = "bin")
}
if (!has_string(sm_status, "GC")) {
  sample_of_interest <- gc_correct(nipt_object = sample_of_interest, method = "bin")
}

n_total <- length(control_group$samples)
n_keep  <- max(10L, round(n_total * 0.8))

control_group <- match_control_group(
  nipt_sample        = sample_of_interest,
  nipt_control_group = control_group,
  mode               = "subset",
  n_of_samples       = n_keep
)

chi_data <- chi_correct(
  nipt_sample        = sample_of_interest,
  nipt_control_group = control_group
)
sample_of_interest <- chi_data$sample
control_group      <- chi_data$control_group

# ---------------------------
# 4) Diagnose control group, remove aberrant controls (robust)
# ---------------------------

control_diag <- diagnose_control_group(nipt_control_group = control_group)

aberrant_names <- character(0)
if (!is.null(control_diag$abberant_scores) && "sample_name" %in% names(control_diag$abberant_scores)) {
  aberrant_names <- unique(control_diag$abberant_scores$sample_name)
}

if (length(aberrant_names)) {
  message("Removing ", length(aberrant_names), " aberrant controls: ", paste(aberrant_names, collapse = ", "))
  for (nm in aberrant_names) {
    control_group <- remove_sample_controlgroup(samplename = nm, nipt_control_group = control_group)
  }
} else {
  message("No aberrant controls detected (or field not present).")
}

# ---------------------------
# 5) Trisomy prediction for ALL requested chromosomes
# ---------------------------

chroms <- parse_chroms(chrom_arg)
chroms <- chroms[chroms >= 1 & chroms <= 22]
if (!length(chroms)) stop("No valid autosomes selected after filtering to 1..22")

n_models     <- 4
n_predictors <- 4
max_elements <- 9

z_results   <- list()
reg_results <- list()

message("Computing Z-score + regression Z-score for chr", min(chroms), "..chr", max(chroms), " ...")

for (chr in chroms) {
  key <- paste0("chr", chr)

  z_results[[key]] <- safe_try(
    calculate_z_score(
      nipt_sample        = sample_of_interest,
      nipt_control_group = control_group,
      chromo_focus       = chr
    )
  )

  reg_results[[key]] <- safe_try(
    perform_regression(
      nipt_sample        = sample_of_interest,
      nipt_control_group = control_group,
      chromo_focus       = chr,
      n_models           = n_models,
      n_predictors       = n_predictors
    )
  )
}

# # NCV across all chromosomes can be VERY slow (brute-force denominator search).
# # Set RUN_NCV_ALL=TRUE if you really want it.
# run_ncv_all <- identical(toupper(Sys.getenv("RUN_NCV_ALL", "FALSE")), "TRUE")
# ncv_templates <- list()
# ncv_results   <- list()

# if (run_ncv_all) {
#   message("RUN_NCV_ALL=TRUE -> Computing NCV templates + scores for all selected chromosomes (may be slow) ...")
#   for (chr in chroms) {
#     key <- paste0("chr", chr)
#     ncv_templates[[key]] <- safe_try(
#       prepare_ncv(
#         nipt_control_group = control_group,
#         chr_focus          = chr,
#         max_elements       = max_elements
#       )
#     )
#     ncv_results[[key]] <- safe_try(
#       calculate_ncv_score(
#         nipt_sample  = sample_of_interest,
#         ncv_template = ncv_templates[[key]]
#       )
#     )
#   }
# } else {
#   message("NCV skipped. To compute NCV for all chromosomes: RUN_NCV_ALL=TRUE Rscript nipter_allchr.R ...")
# }
# ---------------------------
# NCV (default: only 13,18,21; optionally all)
# ---------------------------

run_ncv_all <- identical(toupper(Sys.getenv("RUN_NCV_ALL", "FALSE")), "TRUE")
run_ncv     <- identical(toupper(Sys.getenv("RUN_NCV",     "FALSE")), "TRUE")

ncv_templates <- list()
ncv_results   <- list()

if (run_ncv_all || run_ncv) {
  # which chromosomes to compute NCV for
  ncv_chrom_arg <- Sys.getenv("RUN_NCV_CHRS", "")
  if (run_ncv_all) {
    ncv_chroms <- chroms
  } else if (nzchar(ncv_chrom_arg)) {
    ncv_chroms <- intersect(chroms, parse_chroms(ncv_chrom_arg))
  } else {
    ncv_chroms <- intersect(chroms, c(13, 18, 21))
  }

  message("Computing NCV for: ", paste0("chr", ncv_chroms, collapse = ", "))

  for (chr in ncv_chroms) {
    key <- paste0("chr", chr)

    ncv_templates[[key]] <- safe_try(
      prepare_ncv(
        nipt_control_group = control_group,
        chr_focus          = chr,
        max_elements       = max_elements
      )
    )

    ncv_results[[key]] <- safe_try(
      calculate_ncv_score(
        nipt_sample  = sample_of_interest,
        ncv_template = ncv_templates[[key]]
      )
    )

    # optional debug
    if (inherits(ncv_templates[[key]], "error"))
      message("prepare_ncv failed for ", key, ": ", ncv_templates[[key]]$message)
    if (inherits(ncv_results[[key]], "error"))
      message("calculate_ncv_score failed for ", key, ": ", ncv_results[[key]]$message)
  }
} else {
  message("NCV skipped. Use RUN_NCV=TRUE (default 13/18/21) or RUN_NCV_ALL=TRUE (all selected).")
}
# ---------------------------
# 6) Save results
# ---------------------------

# Full RDS
full_res <- list(
  sample_bam     = sample_bam,
  control_rds    = control_rds,
  chroms         = chroms,
  sample         = sample_of_interest,
  control_group  = control_group,
  z_results      = z_results,
  regression     = reg_results,
  ncv_templates  = ncv_templates,
  ncv_results    = ncv_results,
  diagnostics    = control_diag
)

rds_file <- paste0(out_prefix, ".nipter_allchr_results.rds")
saveRDS(full_res, file = rds_file)

# # Summary TSV (RBZ = mean across models; also output per-model RBZ_1..RBZ_n)
# rbz_sets_mat <- t(sapply(chroms, function(chr) {
#   key <- paste0("chr", chr)
#   zvec <- get_rbz_sets(reg_results[[key]])
#   # pad/truncate to n_models
#   out <- rep(NA_real_, n_models)
#   if (length(zvec)) out[seq_len(min(n_models, length(zvec)))] <- zvec[seq_len(min(n_models, length(zvec)))]
#   out
# }))
# colnames(rbz_sets_mat) <- paste0("RBZ_set", seq_len(n_models))

# summary_df <- data.frame(
#   chr        = chroms,
#   Z_score    = sapply(chroms, function(chr) get_z_score(z_results[[paste0("chr", chr)]])),
#   RBZ_score  = sapply(chroms, function(chr) get_rbz_score(reg_results[[paste0("chr", chr)]], agg = "mean")),
#   rbz_sets_mat,
#   NCV_score  = if (run_ncv_all) sapply(chroms, function(chr) get_ncv_score(ncv_results[[paste0("chr", chr)]])) else NA_real_,
#   stringsAsFactors = FALSE
# )

# tsv_file <- paste0(out_prefix, ".nipter_allchr_summary.tsv")
# write.table(summary_df, file = tsv_file, sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------
# Summary TSV (TRANSPOSED: metrics as rows, chr as columns)
# ---------------------------

# Build per-chromosome labels
chr_labels <- paste0("chr", chroms)  # chr1..chr22 (based on selected chroms)

# Ensure RBZ_set matrix has rownames aligned to chroms
rbz_sets_mat <- t(sapply(chroms, function(chr) {
  key <- paste0("chr", chr)
  zvec <- get_rbz_sets(reg_results[[key]])
  out <- rep(NA_real_, n_models)
  if (length(zvec)) out[seq_len(min(n_models, length(zvec)))] <- zvec[seq_len(min(n_models, length(zvec)))]
  out
}))
colnames(rbz_sets_mat) <- paste0("RBZ_set", seq_len(n_models))

# Original "long" summary (one row per chromosome)
summary_df <- data.frame(
  chr        = chroms,
  Z_score    = sapply(chroms, function(chr) get_z_score(z_results[[paste0("chr", chr)]])),
  RBZ_score  = sapply(chroms, function(chr) get_rbz_score(reg_results[[paste0("chr", chr)]], agg = "mean")),
  rbz_sets_mat,
  NCV_score  = if (run_ncv_all) sapply(chroms, function(chr) get_ncv_score(ncv_results[[paste0("chr", chr)]])) else NA_real_,
  stringsAsFactors = FALSE
)

# Transpose: metrics become rows, chromosomes become columns
mat <- t(as.matrix(summary_df[, setdiff(names(summary_df), "chr"), drop = FALSE]))
colnames(mat) <- chr_labels

# Convert to data.frame with a "metric" column
wide_df <- data.frame(metric = rownames(mat), mat, check.names = FALSE, stringsAsFactors = FALSE)
rownames(wide_df) <- NULL

# Add sample_name row at the top (repeat sample id across chr columns)
sample_id <- basename(out_prefix)
sample_row <- as.data.frame(as.list(c(metric = "sample_name", setNames(rep(sample_id, length(chr_labels)), chr_labels))),
                            check.names = FALSE, stringsAsFactors = FALSE)
wide_df <- rbind(sample_row, wide_df)

# Add chrX and chrY columns (placeholders = NA; sample_name row filled)
wide_df$chrX <- ifelse(wide_df$metric == "sample_name", sample_id, NA)
wide_df$chrY <- ifelse(wide_df$metric == "sample_name", sample_id, NA)

# Reorder columns: metric, chr1..chr22, chrX, chrY
wide_df <- wide_df[, c("metric", chr_labels, "chrX", "chrY"), drop = FALSE]

# Write transposed TSV
tsv_file <- paste0(out_prefix, ".nipter_allchr_summary_wide.tsv")
write.table(wide_df, file = tsv_file, sep = "\t", quote = FALSE, row.names = FALSE)

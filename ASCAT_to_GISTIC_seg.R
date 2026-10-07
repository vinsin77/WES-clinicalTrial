#!/usr/bin/env Rscript
# =============================================================================
# ascat_to_gistic2.R
#
# Converts ASCAT (v3.1.2, WES) outputs into GISTIC2 input files, following
# the method described in [Nature paper, ref 90]:
#   "...the outputs from ASCAT were used to generate input files for GISTIC2,
#    which requires segment coordinates, the number of markers and a
#    log2-scaled copy number. For each segment obtained from segmentation
#    from ASCAT analysis, the normalised depth log ratio was extracted from
#    the ASCAT output files, and the number of assessed loci within each
#    segment was used as the number of markers."
#
# Inputs required:
#   1. segs.rds            - data.frame/data.table with columns:
#                             sample, chr, startpos, endpos, nMajor, nMinor,
#                             nAraw, nBraw  (one row per ASCAT segment)
#   2. Per-sample LogR files - "Tumor.LogR.PCFed.txt"-style files, one per
#                             sample, with rownames formatted as "chr_pos"
#                             and a single column of segmented LogR values.
#                             e.g.:
#                               "1_873343"   0.323207163365363
#                               "1_873344"   0.323207163365363
#
# Outputs:
#   - gistic_segs.txt     : Sample / Chromosome / Start / End / Num_Markers / Seg.CN
#   - gistic_markers.txt  : Marker_ID / Chromosome / Position
#   - conversion_log.txt  : per-sample / per-segment QC summary (markers found,
#                            segments dropped, etc.)
#
# Notes:
#   - Chromosomes are autosomes 1-22 + X (matches the paper's ASCAT run,
#     which included chrX). chrY is excluded (ASCAT/GISTIC2 convention).
#   - Seg.CN here is derived from the RAW segmented LogR value extracted
#     directly from the ASCAT LogR track (NOT recomputed from
#     nMajor/nMinor/nAraw/nBraw - those are ASCAT's integer-rounded calls
#     and are deliberately not used). This matches "normalised depth log
#     ratio...extracted from the ASCAT output files" in the paper.
#   - PLOIDY RESCALING (OPTIONAL - see APPLY_PLOIDY_CORRECTION below):
#     raw ASCAT LogR is expressed relative to the matched normal sample
#     (diploid baseline), not relative to the tumor's own modal ploidy.
#     IMPORTANT: the Nature paper's method ("the normalised depth log ratio
#     was extracted from the ASCAT output files") does NOT mention ploidy
#     rescaling. Applying Seg.CN = LogR - log2(ploidy/2) is a DIFFERENT,
#     additional analytical choice - it assumes the segment with LogR equal
#     to the sample's modal value represents that tumor's own neutral
#     (non-aberrant) state, which is a real biological assumption, not a
#     neutral technical step. This is appropriate if your collaborator
#     specifically wants ploidy-corrected GISTIC2 input, but it means you
#     are running a DIFFERENT analysis than the paper describes, not an
#     exact reproduction of it. Set APPLY_PLOIDY_CORRECTION below to choose:
#         FALSE -> Seg.CN = raw LogR straight from ASCAT (paper reproduction)
#         TRUE  -> Seg.CN = raw LogR - log2(ploidy/2)   (ploidy-corrected)
#     Consider running BOTH versions (toggle + rerun) and comparing results
#     if you're not certain which your collaborator/paper context needs.
#   - Segments with zero overlapping markers in the LogR file are dropped
#     and reported in the log (this usually indicates a coordinate system
#     mismatch - check chr naming / 0- vs 1-based coordinates if it happens
#     a lot).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ----------------------------- CONFIG ---------------------------------------
# >>> EDIT THESE PATHS <<<

segs_rds_path   <- "/~/davidV/AURA/CNV/ASCAT38/segs.RDS"  # combined ASCAT segments file
base_dir        <- "/~/davidV/AURA/CNV/ASCAT38/cnv/"     # contains one subfolder per sample
                                                       # e.g. base_dir/BL1/ascat.Rdata,
                                                       #      base_dir/BL1/Tumor.LogR.PCFed.txt

logr_filename    <- "Tumor.LogR.PCFed.txt"  # filename inside each sample's subfolder
ploidy_filename  <- "ascat.Rdata"           # filename inside each sample's subfolder
                                             # must contain an object called ascat.output
                                             # with a $ploidy field (single named numeric value)

out_dir         <- "/~/myfolder/test_out/cnv_aura/gistic_input"          # where outputs will be written
chroms_keep     <- c(as.character(1:22), "X")  # autosomes + X, matches paper

# >>> ANALYTICAL CHOICE - READ THE NOTE ABOVE BEFORE SETTING THIS <<<
# FALSE = paper-reproduction (raw ASCAT LogR, no ploidy correction)
# TRUE  = ploidy-corrected (Seg.CN = LogR - log2(ploidy/2)) - a different,
#         additional analysis, not what the paper's Methods describes
APPLY_PLOIDY_CORRECTION <- FALSE

# GISTIC2 can run without a custom markers file for WES (it will use the
# segment file's implied marker counts). Set to FALSE to skip generating
# gistic_markers.txt and just omit -mk from your GISTIC2 call.
GENERATE_MARKERS_FILE <- TRUE

# ------------------------------------------------------------------------------

#dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_dir, "conversion_log.txt"), open = "wt")
log_msg <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n")
  writeLines(msg, con = log_con)
}

# ----------------------------- LOAD SEGMENTS ---------------------------------

log_msg("Loading segments from: %s", segs_rds_path)
segs_raw <- readRDS(segs_rds_path)

# segs.RDS is a named list of 113 per-sample data.frames, confirmed structure:
# list(BL1 = data.frame(sample, chr, startpos, endpos, nMajor, nMinor, nAraw, nBraw), ...)
# Flatten into a single data.table. suppressWarnings because rbindlist may warn
# on minor type coercions across samples (character/integer chr column etc.) which
# are harmless - the manual test confirms the output is correct.
if (!is.list(segs_raw) || is.data.frame(segs_raw)) {
  stop("Expected segs.RDS to be a named list of data.frames (one per sample). ",
       "Got: ", class(segs_raw))
}
log_msg("segs.RDS contains %d samples: %s ...",
        length(segs_raw), paste(head(names(segs_raw), 5), collapse = ", "))

segs <- suppressWarnings(rbindlist(segs_raw, use.names = TRUE, fill = TRUE))

required_cols <- c("sample", "chr", "startpos", "endpos")
missing_cols <- setdiff(required_cols, colnames(segs))
if (length(missing_cols) > 0) {
  stop("After flattening segs.RDS, missing required columns: ",
       paste(missing_cols, collapse = ", "),
       "\n  Columns found: ", paste(colnames(segs), collapse = ", "))
}
log_msg("Flattened to %d total segments across %d samples", nrow(segs), length(segs_raw))

# Normalise chr to character without "chr" prefix, for consistent matching
segs[, chr := gsub("^chr", "", as.character(chr))]
segs <- segs[chr %in% chroms_keep]

samples <- sort(unique(segs$sample))
log_msg("Found %d samples, %d segments (after restricting to %s)",
        length(samples), nrow(segs), paste(chroms_keep, collapse = ","))

# ----------------------------- HELPER: PARSE LOGR FILE -----------------------

#' Read one sample's Tumor.LogR.PCFed.txt-style file and return a data.table
#' with columns: chr (character), pos (integer), logr (numeric)
read_logr_file <- function(path) {
  # rownames are like "1_873343", single value column, no guaranteed header
  raw <- fread(path, header = FALSE, sep = "\t")

  # Handle both cases: 2 columns (id, value) or 3 columns (id, value, extra)
  if (ncol(raw) < 2) {
    stop("Unexpected format in LogR file: ", path,
         " (expected at least 2 columns: id, logr)")
  }

  id_col  <- raw[[1]]
  logr_col <- raw[[ncol(raw)]]  # last column = LogR value (robust to optional header col)

  # Strip quotes if present, e.g. "1_873343" -> 1_873343
  id_clean <- gsub('"', "", id_col)

  # Some files may have a header row (e.g. id="" or non-numeric first row) - drop if so
  is_header <- !grepl("^[0-9XYxy]+_[0-9]+$", id_clean)
  if (any(is_header)) {
    n_header <- sum(is_header)
    id_clean <- id_clean[!is_header]
    logr_col <- logr_col[!is_header]
  }

  split_id <- tstrsplit(id_clean, "_", fixed = TRUE)
  chr_vec <- gsub("^chr", "", split_id[[1]])
  pos_vec <- as.integer(split_id[[2]])
  logr_vec <- suppressWarnings(as.numeric(logr_col))

  dt <- data.table(chr = chr_vec, pos = pos_vec, logr = logr_vec)
  dt <- dt[!is.na(pos) & !is.na(logr)]
  setkey(dt, chr, pos)
  dt
}

#' Read a sample's ascat.Rdata and return its continuous ploidy estimate
#' (a single numeric scalar). The object inside is always called
#' "ascat.output" by ASCAT convention, and ascat.output$ploidy is a named
#' numeric vector with one element - we take the value regardless of its
#' name, since ASCAT often names it generically (e.g. "Tumor") rather than
#' with the actual sample ID.
read_sample_ploidy <- function(path) {
  e <- new.env()
  load(path, envir = e)

  if (!exists("ascat.output", envir = e)) {
    stop("No 'ascat.output' object found in: ", path)
  }
  ao <- get("ascat.output", envir = e)

  if (is.null(ao$ploidy) || length(ao$ploidy) == 0) {
    stop("ascat.output$ploidy missing or empty in: ", path)
  }

  as.numeric(ao$ploidy[1])
}

# ----------------------------- MAIN CONVERSION LOOP --------------------------

all_gistic_segs <- vector("list", length(samples))
all_markers     <- vector("list", length(samples))
ploidy_summary  <- vector("list", length(samples))

for (i in seq_along(samples)) {
  s <- samples[i]
  logr_path <- file.path(base_dir, s, logr_filename)

  if (!file.exists(logr_path)) {
    log_msg("[%s] WARNING: LogR file not found at %s - skipping sample", s, logr_path)
    next
  }

  log_msg("[%s] Reading LogR file...", s)
  logr_dt <- read_logr_file(logr_path)
  log_msg("[%s] %d markers loaded", s, nrow(logr_dt))

  # --- Diagnostic: is this file probe-level (many rows per segment, with
  # repeated PCF-constant values) or already one-row-per-segment? Averaging
  # is only meaningful in the former case. Check on the first few segments.
  if (i == 1) {
    check_segs <- segs[sample == s][seq_len(min(5, .N))]
    for (k in seq_len(nrow(check_segs))) {
      cseg <- check_segs[k]
      n_rows <- nrow(logr_dt[chr == cseg$chr & pos >= cseg$startpos & pos <= cseg$endpos])
      log_msg("[%s] DIAGNOSTIC: segment %s:%d-%d has %d LogR rows in file",
              s, cseg$chr, cseg$startpos, cseg$endpos, n_rows)
    }
    log_msg("[%s] DIAGNOSTIC NOTE: if the counts above are mostly 1, the LogR file is",
            s)
    log_msg("  already one-value-per-segment - averaging is a no-op but Num_Markers")
    log_msg("  will be meaningless (always 1) and should NOT be used as-is. If counts")
    log_msg("  are >1 (many probes per segment), the file is probe-level and the")
    log_msg("  current approach (count + mean per segment) is correct. INSPECT THIS")
    log_msg("  before trusting the output.")
  }

  # --- Ploidy rescaling (OPTIONAL - see APPLY_PLOIDY_CORRECTION in CONFIG) ---
  ploidy_offset <- 0  # default: no correction (paper reproduction)
  if (APPLY_PLOIDY_CORRECTION) {
    ploidy_path <- file.path(base_dir, s, ploidy_filename)
    if (!file.exists(ploidy_path)) {
      log_msg("[%s] WARNING: ascat.Rdata not found at %s - skipping sample (cannot rescale without ploidy)",
              s, ploidy_path)
      next
    }
    sample_ploidy <- read_sample_ploidy(ploidy_path)
    ploidy_offset <- log2(sample_ploidy / 2)
    log_msg("[%s] Ploidy = %.3f -> log2(ploidy/2) offset = %.4f applied to all segments",
            s, sample_ploidy, ploidy_offset)
    ploidy_summary[[i]] <- data.table(Sample = s, Ploidy = sample_ploidy, LogR_offset = ploidy_offset)
  } else {
    ploidy_summary[[i]] <- data.table(Sample = s, Ploidy = NA_real_, LogR_offset = 0)
  }

  sample_segs <- segs[sample == s]
  n_total <- nrow(sample_segs)
  n_dropped <- 0L

  seg_rows <- vector("list", n_total)

  for (j in seq_len(n_total)) {
    seg <- sample_segs[j]
    overlapping <- logr_dt[chr == seg$chr & pos >= seg$startpos & pos <= seg$endpos]

    n_markers <- nrow(overlapping)
    if (n_markers == 0) {
      n_dropped <- n_dropped + 1L
      next
    }

    seg_rows[[j]] <- data.table(
      Sample       = s,
      Chromosome   = seg$chr,
      Start        = seg$startpos,
      End          = seg$endpos,
      Num_Markers  = n_markers,
      `Seg.CN`     = mean(overlapping$logr) - ploidy_offset
    )
  }

  seg_rows <- rbindlist(seg_rows, use.names = TRUE)
  all_gistic_segs[[i]] <- seg_rows
  all_markers[[i]] <- logr_dt[, .(chr, pos)]

  log_msg("[%s] %d/%d segments converted, %d dropped (zero overlapping markers)",
          s, nrow(seg_rows), n_total, n_dropped)
}

# ----------------------------- WRITE SEGMENTATION FILE ------------------------

gistic_segs <- rbindlist(all_gistic_segs, use.names = TRUE)
setnames(gistic_segs, "Seg.CN", "Seg.CN")  # keep exact GISTIC2-expected header name

seg_out_path <- file.path(out_dir, "gistic_segs.txt")
fwrite(gistic_segs, seg_out_path, sep = "\t", quote = FALSE)
log_msg("Wrote segmentation file: %s (%d rows, %d samples)",
        seg_out_path, nrow(gistic_segs), length(unique(gistic_segs$Sample)))

# ----------------------------- WRITE MARKERS FILE ------------------------------
# GISTIC2 markers file = union of all marker positions across samples,
# deduplicated, with a unique Marker_ID, Chromosome, Position.
# NOTE: optional for WES - GISTIC2 can run with just -seg and no -mk.
# Set GENERATE_MARKERS_FILE <- FALSE in CONFIG to skip this.

if (GENERATE_MARKERS_FILE) {
  markers_dt <- unique(rbindlist(all_markers, use.names = TRUE))

  # Numeric sort key so chromosome order is 1,2,3...22,X (not lexicographic
  # "1,10,11,...,2,20,..."). chrY is not present since we restricted to
  # chroms_keep earlier, but guard for it anyway.
  markers_dt[, chr_num := fifelse(chr == "X", 23L,
                            fifelse(chr == "Y", 24L,
                                    suppressWarnings(as.integer(chr))))]
  setorder(markers_dt, chr_num, pos)
  markers_dt[, chr_num := NULL]

  markers_dt[, Marker_ID := paste0(chr, "_", pos)]
  markers_out <- markers_dt[, .(Marker_ID, Chromosome = chr, Position = pos)]

  markers_out_path <- file.path(out_dir, "gistic_markers.txt")
  fwrite(markers_out, markers_out_path, sep = "\t", quote = FALSE)
  log_msg("Wrote markers file: %s (%d unique marker positions)",
          markers_out_path, nrow(markers_out))
} else {
  markers_out_path <- NA_character_
  log_msg("GENERATE_MARKERS_FILE = FALSE - skipped markers file. Omit -mk when calling GISTIC2.")
}

# ----------------------------- WRITE PLOIDY QC SUMMARY --------------------------

ploidy_dt <- rbindlist(ploidy_summary, use.names = TRUE)
ploidy_out_path <- file.path(out_dir, "sample_ploidy_offsets.txt")
fwrite(ploidy_dt, ploidy_out_path, sep = "\t", quote = FALSE)
log_msg("Wrote per-sample ploidy/offset QC table: %s", ploidy_out_path)



# Load the segs object
segs <- readRDS("segs.rds")  

# Combine all samples into one data frame
seg_df <- do.call(rbind, segs)

# Build GISTIC 2 seg file (5-column format)
gistic_seg <- data.frame(
  Sample     = seg_df$sample,
  Chromosome = gsub("chr", "", as.character(seg_df$chr)),  # ensure plain "1", "2", ...
  Start      = seg_df$startpos,
  End        = seg_df$endpos,
  Num.Markers = NA,
  Seg.CN     = log2((seg_df$nMajor + seg_df$nMinor) / 2)  # log2 ratio vs diploid
)

# Optional: remove sex chromosomes
gistic_seg <- gistic_seg[!gistic_seg$Chromosome %in% c("X", "Y"), ]

# Optional: cap extreme values (GISTIC can be sensitive to Inf/-Inf)
gistic_seg$Seg.CN[is.infinite(gistic_seg$Seg.CN)] <- NA
gistic_seg <- gistic_seg[!is.na(gistic_seg$Seg.CN), ]

# Write output
write.table(gistic_seg, file = "all_samples_gistic.seg",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done:", nrow(gistic_seg), "segments,", length(unique(gistic_seg$Sample)), "samples\n")

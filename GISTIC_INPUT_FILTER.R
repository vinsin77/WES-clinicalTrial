#ASCAT was run on the whole cohort ( Baseline and  Surgery together), but 5 samples were not removed
#Before GISTIC2, I excluded 5 (here is 4) samples that were already removed from the metadata
#And I will run the GISTIC2 for baseline samples only


module load R/4.4.2-gfbf-2024a
library(data.table)
in_dir  <- "gistic_input"
out_dir <- "gistic_input_BL"
exclude_ids <- c("BL1", "BL31", "BL86", "BL88")
segs <- fread(file.path(in_dir, "gistic_segs.txt"))
# typo check: every excluded ID must exist
stopifnot(all(exclude_ids %in% segs$Sample))

keep BL samples, then drop the excluded ones
segs_bl <- segs[startsWith(Sample, "BL") & !Sample %in% exclude_ids]
cat("Samples before:", uniqueN(segs$Sample),
    "| after:", uniqueN(segs_bl$Sample), "\n")

fwrite(segs_bl, file.path(out_dir, "gistic_segs_BL.txt"), sep = "\t", quote = FALSE)
fwrite(data.table(Sample = sort(unique(segs_bl$Sample))),
       file.path(out_dir, "samples_used.txt"), quote = FALSE)

# markers file is unchanged (same array), so just copy it
file.copy(file.path(in_dir, "gistic_markers.txt"),
          file.path(out_dir, "gistic_markers.txt"))

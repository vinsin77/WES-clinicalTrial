
#
# GISTIC command; then 
# gistic_output/all_lesions.conf_95.txt; then
# R script (post-processing) #this is where we apply the  thresholds

#working with absolute copy number values

# Deep deletion    : CN ≤ 0.7
# Loss             : 0.7 < CN ≤ 1.6
# Neutral          : 1.6 < CN ≤ 2.5
# Gain             : 2.5 < CN ≤ 5
# Amplification    : CN > 5

library(dplyr)

# 
# Read GISTIC all_lesions file
all_lesions <- read.delim("gistic_output/all_lesions.conf_95.txt", 
                           header = TRUE, check.names = FALSE)



# Filter CN values rows only

cn_rows <- all_lesions[all_lesions$`Amplitude Threshold` == "Actual Copy Change Given", ]


# Build CN matrix and normalize to ploidy 2

cn_matrix <- apply(cn_rows[, 10:ncol(cn_rows)], 2, as.numeric) + 2
cat("CN matrix dimensions:", dim(cn_matrix), "\n")

# Classify using paper's thresholds
classify_cn <- function(x) {
  case_when(
    x <= 0.7 ~ "Deep deletion",
    x <= 1.6 ~ "Loss",
    x <= 2.5 ~ "Neutral",
    x <= 5.0 ~ "Gain",
    x >  5.0 ~ "Amplification",
    TRUE     ~ NA_character_
  )
}

status_matrix <- apply(cn_matrix, c(1,2), classify_cn)
rownames(status_matrix) <- cn_rows$Descriptor

# Quick summary:how many peaks per category

print(table(status_matrix))
#Amplification          Gain          Loss       Neutral 
          #144           973           978          2952


#  Save outputs

# Status matrix (Deep deletion / Loss / Neutral / Gain / Amplification)
write.csv(status_matrix, "cn_status_classified.csv")

# CN matrix (absolute values, ploidy-normalized)
cn_matrix_df <- as.data.frame(cn_matrix)
rownames(cn_matrix_df) <- cn_rows$Descriptor
write.csv(cn_matrix_df, "cn_absolute_values.csv")


cat("  cn_absolute_values.csv    — absolute CN values per peak per sample\n")


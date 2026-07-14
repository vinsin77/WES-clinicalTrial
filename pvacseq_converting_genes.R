#this step is required for pvaqseq (neoantigen calculation); before adding the expression data to VCF
##The TPM file uses gene symbols (A1BG, A1BG-AS1, 7SK etc.) as row names, but pvacseq/VEP expects Ensembl gene IDs (ENSG00000...). 
##But vcf-expression-annotator needs the IDs to match the Gene field in your VEP CSQ annotation.
## We need to convert gene symbols to Ensembl IDs first:

library(biomaRt)

# Load and convert to data frame
data_tpm <- readRDS("RNAtpm.RDS") #expression file
data1 <- as.data.frame(data_tpm)

# Extract BL-100 only
bl100_tpm <- data.frame(
  gene_symbol = rownames(data1),
  TPM = data1[, "BL-100"]
)

head(bl100_tpm)
nrow(bl100_tpm)

# Connect to Ensembl and map symbols to Ensembl IDs
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

mapping <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = bl100_tpm$gene_symbol,
  mart = mart
)

# Merge
bl100_mapped <- merge(bl100_tpm, mapping, by.x = "gene_symbol", by.y = "hgnc_symbol")

# Keep only Ensembl ID and TPM columns
bl100_final <- bl100_mapped[, c("ensembl_gene_id", "TPM")]
colnames(bl100_final) <- c("gene_id", "TPM")

# Remove duplicates if any
bl100_final <- bl100_final[!duplicated(bl100_final$gene_id), ]

# Write to file
write.table(bl100_final, "BL100_tpm.tsv", sep="\t", row.names=FALSE, quote=FALSE)

head(bl100_final)
nrow(bl100_final)

#!! Now annotate the VCF with the expression data !!

library(biomaRt)
data_tpm <- readRDS("RNAtpm.RDS")
data1 <- as.data.frame(data_tpm)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id"),
    filters = "hgnc_symbol",
    values = rownames(data1),
    mart = mart
)
for (sample in colnames(data1)) {

    sample_tpm <- data.frame(
        gene_symbol = rownames(data1),
        TPM = data1[, sample]
    )

    sample_mapped <- merge(
        sample_tpm,
        mapping,
        by.x = "gene_symbol",
        by.y = "hgnc_symbol"
    )

    sample_final <- sample_mapped[, c("ensembl_gene_id", "TPM")]
    colnames(sample_final) <- c("gene_id", "TPM")

    sample_final <- sample_final[!duplicated(sample_final$gene_id), ]

    write.table(
        sample_final,
        file = file.path("TPM_files", paste0(gsub("-", "", sample), "_tpm.tsv")),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
}

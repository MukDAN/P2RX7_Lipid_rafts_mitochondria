setwd("/media/mukesh/HIKVISION/mbcd_rna")
countData <- read.table("counts/all_samples_counts.txt",
                        header=TRUE,
                        row.names=1,
                        comment.char="#")
# Keep only count columns
countData <- countData[,6:ncol(countData)]

colnames(countData)
coldata <- data.frame(
  row.names = colnames(countData),
  condition = c("control","control","control",
                "treated","treated","treated")
)
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = coldata,
                              design = ~ condition)

dds <- DESeq(dds)

res <- results(dds)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Convert DESeq2 result to dataframe
res_df <- as.data.frame(res)

# Add ENSEMBL IDs
res_df$ENSEMBL_ID <- rownames(res_df)

# Remove version numbers
gene_ids <- sub(
  "\\..*",
  "",
  res_df$ENSEMBL_ID
)

# Add gene symbols
res_df$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Reorder columns
res_df <- res_df[
  ,
  c(
    "ENSEMBL_ID",
    "SYMBOL",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj"
  )
]

# Save all genes
write.csv(
  res_df,
  "DESeq2_all_genes.csv",
  row.names = FALSE
)


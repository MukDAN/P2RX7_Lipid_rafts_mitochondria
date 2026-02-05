library(clusterProfiler)
library(org.Hs.eg.db)   # Human gene annotation
# ===============================
# 2. Load data
# ===============================
df1 <- read.csv("Lipid_raft_associated_proteins - Sheet1.csv", stringsAsFactors = FALSE)%>%as.data.frame()
df2 <- read.csv("raft_data_JJP_proteomic.csv", stringsAsFactors = FALSE)%>%as.data.frame()
df3 <- read.csv("raft_data_proteomic.csv", stringsAsFactors = FALSE)%>%as.data.frame()

colnames(df1)

# ===============================
# 3. Harmonize column names
# ===============================
df1 <- df1 %>%
  dplyr::rename(UniProt = Protein_Accession,
                GeneName = Gene_Name)

df2 <- df2 %>%
  dplyr:: select(-X) %>%
  dplyr::rename(UniProt = UNIPROT,
                GeneName = Gene.name)

df3 <- df3 %>%
  dplyr::select(-X) %>%
  dplyr::rename(UniProt = UNIPROT,
                GeneName = Gene.Symbol)

# ===============================
# 4. Merge + clean
# ===============================
all_proteins <- bind_rows(df1, df2, df3) %>%
  filter(!is.na(GeneName), GeneName != "") %>%
  mutate(GeneName = str_trim(GeneName))

# Remove duplicates (final protein list)
all_proteins_unique <- all_proteins %>%
  distinct(GeneName, .keep_all = TRUE)

# Use unique UniProt IDs
uniprot_ids <- all_proteins_unique$UniProt
uniprot_ids <- uniprot_ids[!is.na(uniprot_ids) & uniprot_ids != ""]
length(uniprot_ids)



library(clusterProfiler)
library(org.Hs.eg.db)
library(clusterProfiler)
library(org.Hs.eg.db)

# UniProt → Entrez
gene_map <- bitr(
  uniprot_ids,
  fromType = "UNIPROT",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

ego_bp <- enrichGO(
  gene          = gene_map$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_cc <- enrichGO(
  gene          = gene_map$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego_cc)

ego_mf <- enrichGO(
  gene          = gene_map$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

# KEGG enrichment
ego_kegg <- enrichKEGG(
  gene         = gene_map$ENTREZID,
  organism     = "hsa",          # 'hsa' for human, 'mmu' for mouse
  pvalueCutoff = 0.05,
  pAdjustMethod= "BH"
)

# Optional: make gene symbols readable
#ego_kegg <- setReadable(ego_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# View results


head(ego_mf)

bp_df    <- as.data.frame(ego_bp)
mf_df    <- as.data.frame(ego_mf)
cc_df    <- as.data.frame(ego_cc)
kegg_df  <- as.data.frame(ego_kegg)
bp_df$logAdjP   <- -log10(bp_df$p.adjust)
mf_df$logAdjP   <- -log10(mf_df$p.adjust)
cc_df$logAdjP   <- -log10(cc_df$p.adjust)
kegg_df$logAdjP <- -log10(kegg_df$p.adjust)

library(ggplot2)

plot_enrich <- function(df, ylab, title) {
  ggplot(
    df,
    aes(
      x = logAdjP,
      y = reorder(Description, logAdjP),
      fill = logAdjP
    )
  ) +
    geom_col() +
    scale_fill_viridis_c(option = "plasma") +
    labs(
      x = expression(-log[10]("Adjusted p-value")),
      y = ylab,
      title = title,
      fill = "logAdjP"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 9)
    )
}

library(dplyr)

top_bp <- bp_df %>%
  slice_max(order_by = logAdjP, n = 20)

top_mf <- mf_df %>%
  slice_max(order_by = logAdjP, n = 20)

top_cc <- cc_df %>%
  slice_max(order_by = logAdjP, n = 20)

top_kegg <- kegg_df %>%
  slice_max(order_by = logAdjP, n = 20)

p_bp <- plot_enrich(
  top_bp,
  "GO Biological Process",
  "Top 20 enriched GO-BP"
)

p_mf <- plot_enrich(
  top_mf,
  "GO molecular function",
  "Top 20 enriched-MF"
)


p_cc <- plot_enrich(
  top_cc,
  "GO cellular components",
  "Top 20 enriched GO-CC "
)

p_kegg <- plot_enrich(
  top_kegg,
  "KEGG pathways",
  "Top 20 enriched KEGG pathways"
)

p_kegg
library(dplyr)


library(patchwork)

(p_bp | p_mf) /
  (p_cc | p_kegg) &
  theme(
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x  = element_text(size = 16,face="bold"),
    axis.text.y  = element_text(size = 16, face="bold"),
    plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 16)
  )

# Add a column indicating the ontology/source
bp_df <- bp_df %>% mutate(Category = "BP")
mf_df <- mf_df %>% mutate(Category = "MF")
cc_df <- cc_df %>% mutate(Category = "CC")
kegg_df <- kegg_df %>% mutate(Category = "KEGG")

# Combine all into a single data frame
combined_df1 <- bind_rows(bp_df, mf_df, cc_df, kegg_df)

# View as table
combined_df
write.csv(combined_df1,"Top20_GO_KEGG_analysis_of Rafts_proteins.csv")


top_bp <- bp_df %>%
  slice_max(order_by = logAdjP, n = 20)%>%mutate(Category = "BP")

top_mf <- mf_df %>%
  slice_max(order_by = logAdjP, n = 20)%>%mutate(Category = "MF")

top_cc <- cc_df %>%
  slice_max(order_by = logAdjP, n = 20)%>%mutate(Category = "CC")

top_kegg <- kegg_df %>%
  slice_max(order_by = logAdjP, n = 20)%>%mutate(Category = "KEGG")
write.csv(combined_df,"GO_KEGG_analysis_of Rafts_proteins.csv")

setwd("/media/mukesh/84EABC40EABC2FF2/Paper_2/")
# Load required libraries
library(readxl)
library(tidyverse)
library(STRINGdb)
library(igraph)
library(RCy3)
library(biomaRt)
library(jsonlite)
library(stringr)
library(STRINGdb)
library(igraph)
library(ggraph)
library(ggplot2)

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
library(readxl)
mitocarta<-read_xlsx("/media/mukesh/84EABC40EABC2FF2/Paper_2/main_figure/final/main_figure/revision/supllimentary/mitocarta.xlsx")

all_protein<-all_proteins_unique%>%filter(UniProt%in%c(mitocarta$UniProt))
mitocarta_lipid_rafts<-mitocarta%>%filter(UniProt%in%c(all_protein$UniProt))

str(mitocarta)
library(dplyr)
library(tidyr)
library(stringr)

mito1 <- mitocarta %>%
  separate_rows(MitoCarta3.0_MitoPathways,
                sep=" \\| ") %>%
  mutate(
    Major_pathway=
      str_trim(
        str_extract(
          MitoCarta3.0_MitoPathways,
          "^[^>]+"
        )
      )
  )

mito1_gene_unique <- mito1[!duplicated(mito1$Symbol), ]
library(dplyr)

mito1_gene_unique <- mito1 %>%
  distinct(Symbol, .keep_all = TRUE)


length(unique(mito1_gene_unique))
summary_tbl <- mito1 %>%
  count(Major_pathway,
        sort=TRUE)

summary_tbl
summary_tbl2 <- mito1 %>%
  filter(!Major_pathway %in% c("0"),
         !is.na(Major_pathway)) %>%
  count(Major_pathway, sort = TRUE)

summary_tbl2 %>%
  mutate(
    Percent = round(
      n/sum(n)*100,
      1
    )
  )

library(ggplot2)

ggplot(
  summary_tbl2,
  aes(
    x = reorder(Major_pathway,n),
    y = n
  )
)+
  geom_col(width=.8)+
  coord_flip()+
  theme_bw(base_size=14)+
  labs(
    x="",
    y="Number of proteins",
    title="Functional distribution of lipid raft-associated mitochondrial proteins"
  )

oxphos <- mito1 %>%
  filter(
    Major_pathway=="OXPHOS"
  )


table<-
  
  table(
  sub(
    ".*Complex ([IVX]+).*",
    "Complex \\1",
    oxphos$MitoCarta3.0_MitoPathways
  )
)

library(dplyr)

summary_tbl2 <- summary_tbl2 %>%
  mutate(
    Percent = round(n/sum(n)*100, 1),
    label = paste0(Major_pathway,
                   "\n",
                   Percent,
                   "%")
  )

summary_tbl2
library(ggplot2)
summary_tbl2 <- summary_tbl2 %>%
  filter(
    !is.na(Major_pathway),
    Major_pathway != "0"
  )

ggplot(summary_tbl2,
       aes(x = 2,
           y = n,
           fill = Major_pathway)) +
  
  geom_col(width = 1,
           color = "white") +
  
  coord_polar(theta = "y") +
  
  xlim(0.5, 2.5) +
  
  geom_text(
    aes(label = paste0(Percent, "%")),
    position = position_stack(vjust = 0.5),
    size = 8
  ) +
  
  theme_void() +
  
  labs(
    title = " ",
    fill = "Major pathway"
  ) +
  
  theme(
    legend.title = element_text(
      size = 18,
      face = "bold"
    ),
    legend.text = element_text(
      size = 16,
      face = "bold"
    )
  )
  

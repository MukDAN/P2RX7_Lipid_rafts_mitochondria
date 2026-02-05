setwd("/media/mukesh/84EABC40EABC2FF2/Paper_2/")
library(igraph)
library(data.table)
mapped<-read.csv("lipid_raft_data_for_gene_entich.csv")
# ===============================
# 6. Retrieve and clean PPI network
# ===============================
ppi_all <- fread("9606.protein.links.v11.5.txt.gz")
valid_ids <- unique(mapped$STRING_id)

ppi_edges <- ppi_all[
  protein1 %in% valid_ids &
    protein2 %in% valid_ids
]


# Confidence filter (recommended)
ppi_edges <- ppi_edges[combined_score >= 400]
ppi_edges%>%filter(protein1=="9606.ENSP00000330696")
ppi_edges%>%filter(protein2=="9606.ENSP00000330696")  #P2RX7

# Safety check
ppi_edges <- ppi_edges[!is.na(protein1) & !is.na(protein2)]

# Sanity check
stopifnot(nrow(ppi_edges) > 0)

# Build igraph network


#ppi_edges <- string_db$get_interactions(mapped$STRING_id)

ppi_edges_clean <- ppi_edges %>%
  distinct(protein1, protein2, combined_score)


# ===============================
# 7. Build igraph object
# ===============================
g <- graph_from_data_frame(
  d = ppi_edges_clean,
  directed = FALSE
)


# Number of nodes (proteins)
vcount(g)

# Number of edges (PPIs)
ecount(g)

# Load required libraries
library(igraph)
library(ggraph)
library(ggplot2)
library(tidyverse)

# -------------------------
# Compute network properties
# -------------------------
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g, normalized = TRUE)
V(g)$closeness <- closeness(g, normalized = TRUE)
V(g)$eigenvector <- eigen_centrality(g)$vector
V(g)$cluster_coeff <- transitivity(g, type = "local")

# Shortest paths summary
avg_path_length <- average.path.length(g)
cat("Average shortest path length:", avg_path_length, "\n")

# -------------------------
# ===============================
# Libraries
# ===============================
library(igraph)
library(dplyr)
library(tidyr)
library(ggplot2)

# ===============================
# 1. Extract node-level metrics
# ===============================
node_data <- data.frame(
  name          = V(g)$name,
  degree        = V(g)$degree,
  betweenness   = V(g)$betweenness,
  closeness     = V(g)$closeness,
  eigenvector   = V(g)$eigenvector,
  cluster_coeff = V(g)$cluster_coeff
)

# ===============================
# 2. Convert to long format
# ===============================
node_data_long <- node_data %>%
  pivot_longer(
    cols      = -name,
    names_to  = "Metric",
    values_to = "Value"
  )

# ===============================
# 3. Min–max scaling per metric
# ===============================
node_data_scaled <- node_data_long %>%
  group_by(Metric) %>%
  mutate(
    Value_scaled = (Value - min(Value, na.rm = TRUE)) /
      (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE))
  ) %>%
  ungroup()

# ===============================
# 4. Dark scientific color palette
# ===============================
dark_palette <- c(
  degree        = "#1B4F72",  # deep blue
  betweenness   = "#7D3C98",  # dark purple
  closeness     = "#117864",  # dark teal
  eigenvector   = "#B03A2E",  # dark red
  cluster_coeff = "#7E5109"   # dark brown
)

# ===============================
# 5. Plot: Dark violin + jitter
# ===============================
ggplot(
  node_data_scaled,
  aes(x = Metric, y = Value_scaled, fill = Metric, color = Metric)
) +
  
  # Violin
  geom_violin(
    alpha     = 0.5,
    color     = "black",
    linewidth = 0.6,
    trim      = TRUE
  ) +
  
  # Jitter points
  geom_jitter(
    width = 0.2,
    size  = 3,
    alpha = 0.9
  ) +
  
  # Colors
  scale_fill_manual(values = dark_palette) +
  scale_color_manual(values = dark_palette) +
  
  # Labels
  labs(
    title    = "a",
    y = "Scaled value (0–1)",
    x = ""
  ) +
  
  # Theme (publication-ready)
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size  = 30,
      face  = "bold"
    ),
    
    axis.text.y = element_text(
      size = 30,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 30,
      face = "bold"
    ),
    
    plot.title = element_text(
      size = 30,
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 18,
      hjust = 0.5
    ),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme(
    plot.title = element_text(
      size = 46,
      face = "bold",
      hjust = 0   # ← extreme left
    ),
    
    plot.subtitle = element_text(
      size = 18,
      hjust = 0
    )
  )

#############
library(igraph)

hub_ids <- sub("^9606\\.", "", names(hub_nodes))

library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

hub_map <- getBM(
  attributes = c("ensembl_peptide_id", "hgnc_symbol"),
  filters = "ensembl_peptide_id",
  values = hub_ids,
  mart = mart
)

hub_df <- data.frame(
  ensembl_peptide_id = hub_ids,
  Degree = as.numeric(hub_nodes)
)
# Assuming 'g' is your PPI network
bet <- betweenness(g, normalized = TRUE)
V(g)$size <- bet*5  # scale node size by betweenness
V(g)$color <- ifelse(bet > quantile(bet,0.9), "red", "skyblue")  # top 10% as bottlenecks

plot(g, vertex.label=V(g)$name, vertex.size=V(g)$size, vertex.color=V(g)$color)



#### 
# Gene ontology enrichment
# Author: Filippa Bertilsson
# Contact: bertilsson.filippa@igp.uu.se
####

library(stringr)
# Load packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("fgsea")
library(fgsea)
BiocManager::install("clusterProfiler")
library("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
library(tidyverse)
library(ggplot2)
library(viridis)
BiocManager::install("GO.db", force = TRUE)
library(GO.db)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(forcats)
library(readxl)
library(writexl)
set.seed(42)
select <- dplyr::select
slice <- dplyr::slice
rename <- dplyr::rename
filter <- dplyr::filter


setwd("~/Library/CloudStorage/Box-Box/Cilia/Publication/Data/Data")


# Gene ontology

# Load in clusters from the heatmap
ALL.cluster.GO.7 <- readRDS("ALL Cluster GO 7 CLUSTERS 2025-02-13 MEDIAN.rds") 

# Prep lists for analysis
ego_results_ALL <- list()
clusters <- unique(ALL.cluster.GO.7$cluster)


# Loop to identify GO terms
for (cluster_name in clusters) {
  cluster_data <- ALL.cluster.GO.7 %>% filter(cluster == cluster_name)
  
  # Settings
  ego_results_ALL[[cluster_name]] <- enrichGO(
    gene           = cluster_data$Protein, # Proteins we look at
    OrgDb          = org.Hs.eg.db, # Organism data base (Human)
    keyType        = "SYMBOL", # Specifies the type of ID (gene symbol)
    ont            = "ALL", # Ontology: All, includes BP, MF, and CC
    pAdjustMethod  = "BH", # Adjustment method for multiple testing (BH = Benjamini-Hochberg)
    pvalueCutoff   = 0.1, # Keeps only GO terms with unadjusted p-value ≤ 0.1
    qvalueCutoff   = 0.2 # Keeps only GO terms with adjusted p-value ≤ 0.2
  )
  
  ego_results_ALL[[cluster_name]] <- clusterProfiler::simplify(
    ego_results_ALL[[cluster_name]], 
    cutoff = 0.7, # Terms with ≥70% overlap in gene sets are considered redundant
    by = "p.adjust", # Among redundant terms, we keep the one with the lowest adjusted p-value
    select_fun = min # How to choose among redundant terms (here, the one with minimum p.adjust is chosen)
  )
}


# Make list to vector
combined_data_ALL <- vector("list", length(ego_results_ALL))

# Organize GO results
for (cluster_name in names(ego_results_ALL)) {  # Use annotated cluster names
  result <- ego_results_ALL[[cluster_name]]
  if (nrow(result) > 0) {
    combined_data_ALL[[cluster_name]] <- data.frame(
      Cluster = rep(cluster_name, length(result$Description)),  # Annotated cluster name
      Term = result$Description,
      Gene = result$geneID,
      Ontology = result$ONTOLOGY,
      ID = result$ID,
      Ratio = result$GeneRatio,
      Count = result$Count,
      p.adjust = result$p.adjust
    )
  } 
}

ALL.GO.DATA <-  do.call(rbind, combined_data_ALL)



# Filter out terms with p-value < 0.05
ALL.GO.DATA.f <- ALL.GO.DATA %>% dplyr::filter(p.adjust < 0.05)

# Calculate rich factor
ALL.GO.DATA.f <- ALL.GO.DATA %>%
  mutate(richFactor = Count / as.numeric(sub(".*/", "", Ratio)))

# Put in cluster names from heatmap
clusters_order_ALL <- c("Cluster.7", "Cluster.1", "Cluster.5",  "Cluster.3", "Cluster.2", "Cluster.6", "Cluster.4")

# Choose only the top 3 terms for each cluster to then display in the plot
data.ALL.count <- ALL.GO.DATA.f %>%
  group_by(Cluster) %>%
  top_n(n = 3, wt = Count) %>%
  arrange(Cluster, p.adjust) %>%
  group_by(Cluster) %>%
  slice_head(n = 3) %>%
  ungroup()

data.ALL.count <- data.ALL.count %>%
  filter(Cluster %in% clusters_order_ALL) %>%
  mutate(Cluster = factor(Cluster, levels = clusters_order_ALL))

# Add line breaks to GO terms (split into ~20-character chunks)
data.ALL.count$Term_wrapped <- str_wrap(data.ALL.count$Term, width = 30)


# Lollipop plot
# Settings for the plot
p.ALL.COUNT <- ggplot(data.ALL.count, aes(y = fct_reorder(Term_wrapped, richFactor), x = richFactor)) +
  geom_segment(aes(y = Term_wrapped, yend = Term_wrapped, x = 0, xend = richFactor)) +  # Horizontal lines
  geom_point(aes(color = p.adjust, size = Count)) +  # Dots at the end of the lines
  scale_color_gradientn(
    colours = c("#FADADD", "#FFA07A", "#D64550"),
    trans = "log10",
    guide = guide_colorbar(reverse = TRUE, order = 1)
  ) +
  scale_size_continuous(range = c(2, 10)) +
  facet_wrap(~ Cluster, ncol = 1, scales = "free_y", strip.position = "right") +  
  labs(y = " ", x = "Rich Factor") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12), 
    strip.background = element_blank(),
    strip.text = element_text(size = 14, angle = -90, hjust = 0.5),  
    strip.placement = "outside",  
    panel.spacing = unit(1, "lines"),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()   
  )

# Print the plot
print(p.ALL.COUNT) # save as B_LollipopGO.pdf (7.39x6.58)  

# FYI 
# Clusters in the plot are now displayed in this order from top to bottom: "Cluster.7", "Cluster.1", "Cluster.5",  "Cluster.3", "Cluster.2", "Cluster.6", "Cluster.4"
# But for the figures in the paper we renamed them to be in a logical order form top to bottom like this:
# Original cluster => Our cluster name
# Cluster 7 => Cluster 1
# Cluster 1 => Cluster 2
# Cluster 5 => Cluster 3
# Cluster 3 => Cluster 4
# Cluster 2 => Cluster 5
# Cluster 6 => Cluster 6
# Cluster 4 => Cluster 7








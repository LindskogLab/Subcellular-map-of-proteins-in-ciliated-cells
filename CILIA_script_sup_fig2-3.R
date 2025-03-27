
#### 
# Supplimentary figures 2 and 3
# Author: Loren Méar
####

# load packages 
library(tidyverse)
library(data.table)
library(readxl)
library(pheatmap)

# REMEMBER
# Clusters in the plots are now displayed in this order from top to bottom: "Cluster.7", "Cluster.1", "Cluster.5",  "Cluster.3", "Cluster.2", "Cluster.6", "Cluster.4"
# But for the figures in the paper we renamed them to be in a logical order form top to bottom like this:
# Original cluster => Our cluster name
# Cluster 7 => Cluster 1
# Cluster 1 => Cluster 2
# Cluster 5 => Cluster 3
# Cluster 3 => Cluster 4
# Cluster 2 => Cluster 5
# Cluster 6 => Cluster 6
# Cluster 4 => Cluster 7

# Aesthetic 
ann_colors <- c("1" = "#F4ACB7", "2" = "#C0DFA1", "3" = "#C19EE0",
                "4" = "#FCBF6E", "5" = "#F08080", "6" = "#FFF3B0", "7" = "#90C2E7")

fe_colors <- c( "1" = "#c7c9c4", "2" = "#8f918f",
                "3" = "#946b91", "4" = "#50364D", "5" = "#241826" )

################################################################################
#### Supplemental figure 2
################################################################################

setwd()

# list of the 187 candidates
ListProt <- (combined_matrix_table)

# scRNA data from the atlas _ 13/02/2025
scRNA_Atlas <- fread() # scRNA data from HPA (proteinatlas.org)
Cilia_scRNA_Atlas <- left_join(ListProt, scRNA_Atlas, by = join_by("Gene"=="Gene name"  ))

# selected cell type to include 
celltype <- read_excel() # File with the cell types we wanted to display
selected_cell_types <- celltype %>% 
  filter(Keep == "yes") %>% 
  pull(`Cell Type`)

# Data to plot
heatmap_matrix <- Cilia_scRNA_Atlas %>%
  filter(`Cell type` %in% selected_cell_types) %>% 
  dplyr::select(Gene, `Cell type`, nTPM)%>%
  mutate(nTPM_log = log10(nTPM + 1)) %>%
  dplyr::select(Gene, `Cell type`, nTPM_log) %>%
  pivot_wider(names_from = `Cell type`, values_from = nTPM_log, values_fill = 0) %>%
  column_to_rownames(var = "Gene") %>%
  as.matrix()


# heatmap
pheatmap(heatmap_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "none",  
         show_rownames = F, 
         show_colnames = TRUE, 
         color = colorRampPalette(c("white", "#FFD3B6", "#FF9F80", "#FF6F61", "#B03040"))(100), 
         fontsize_col = 10, 
         angle_col = 45)



################################################################################
#### Supplemental figure 3
################################################################################
setwd()

# Data from HPP portal downloaded the 13/02/2025 : score FE and PE for the 187 candidates 
score_187 <- fread()

# Data overview 
table(score_187$`PE score`, score_187$`FE score`)

# Cluster information 
cluster.info <- read_excel() %>% select("Protein", "cluster")

# Merge data (HPP_Cluster)
score_187_cl <- score_187 %>%
  left_join(cluster.info, by = join_by("Gene" == "Protein")) %>%
  mutate(
    `PE score` = as.factor(`PE score`),
    `FE score` = as.factor(`FE score`),
    cluster = as.factor(cluster)
  ) %>%
  rename_with(~ gsub(" ", "_", .x))

# Data overview
table(score_187_cl$FE_score, score_187_cl$cluster)

# Data frame to plot 
df_plot <- score_187_cl %>%
  dplyr::count(cluster, FE_score) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) 


# Plot setting
p <- ggplot(df_plot, aes(x = factor(cluster), y = prop, fill = factor(FE_score))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "white", width = 0.7) +
  scale_fill_manual(values = fe_colors, name = "FE score") +
  theme_minimal() +
  labs(x = "", y = "Proportion") +
  facet_wrap2(~ cluster, scales = "free_x", nrow = 1, strip = strip_themed(
    background_x = elem_list_rect(fill = unname(ann_colors[as.character(unique(df_plot$cluster))]))
  )) +
  theme(
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

# Print plot
print(p)

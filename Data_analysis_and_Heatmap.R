
#### 
# Data analysis of pixel values to determine subspatial localization of proteins in ciliated cells
# Author: Filippa Bertilsson
# Contact: bertilsson.filippa@igp.uu.se
####



# Load packages
library(readr)
library(stringr)
library(tools)
library(writexl)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(readxl)
library(pheatmap)
library(purrr)
theme_set(theme_minimal())

# Specify packages
select <- dplyr::select
slice <- dplyr::slice
rename <- dplyr::rename
filter <- dplyr::filter


# Set working directory
setwd("/Your/Directory/Here") 
dir()


### Find the SlideIDs and put them in a list
# Specify the path to your folder containing the CSV files
folder_path <- "/Your/Folder/Path/Here" 

# Get the list of file names in the folder
# Read in whole file names
file_names <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Extract SlideIDs from file names
slide_ids <- unique(sapply(file_names, function(file) {
  basename(file) %>%
    str_extract(pattern = "(?<=_)[A-Za-z0-9]+(?=\\.csv$)")
}))

# Define the list of tissues and SlideIDs
tissues <- c("Bronchus", "Nasopharynx", "Cervix", "Endometrium", "Fallopian") # The tissues we used in this study
SlideIDs <- c(slide_ids)

# Make SlideID list into table
gene_table <- as_data_frame(SlideIDs) 

# Define core names for each tissue
core_names <- list(
  Bronchus = c("7C", "7D", "8C", "8D", "9C", "9D"),
  Nasopharynx = c("4C", "4D", "5C", "5D", "6C", "6D"),
  Cervix = c("1C", "1D", "2C", "2D", "3C", "3D"),
  Endometrium = c("4A", "4B", "5A", "5B", "6A", "6B"),
  Fallopian = c("1A", "1B", "2A", "2B", "3A", "3B")
)


# Read in meta data
# Here, we read in information about which images (core names) were removed during manual QC (from our internal LIMS)
# This information is later used in the loop where we process the data to exclude results from images that didn't pass QC
MetaData.Cilia <- read_excel("file_here") %>% mutate(`Slide ID` = as.character(`Slide ID`))

max_splits_Cilia <- max(sapply(strsplit(MetaData.Cilia$discarded_images, ","), length))

new_cols_Cilia <- lapply(strsplit(MetaData.Cilia$discarded_images, ","),
                         function(x) c(x, rep(NA, max_splits_Cilia - length(x))))

cores_to_delete <- cbind(MetaData.Cilia, do.call(rbind, new_cols_Cilia)) %>%
  filter(State != "Failed", State != "Suggested")  # This removes proteins in the staining stages Failed and Suggested (internal system)

colnames(cores_to_delete)[(ncol(cores_to_delete) - max_splits_Cilia + 1):ncol(cores_to_delete)] <-
  paste("cores_to_delete", 1:max_splits_Cilia, sep = "_")
# Now there is a column for each core that should be removed when we read in the data


# Data processing
# Create an empty list to store the results
result_list_pixel <- list()

# Ensure that the cores_to_delete file has a column for gene names (assumed "Gene")
# Check this part to really see what it does
if (!"Gene" %in% colnames(cores_to_delete)) {
  stop("The cores_to_delete file must have a 'Gene' column.")
}

# Loop through results form image analysis to collect everyting in one list
for (tissue in tissues) {
  for (slide_id in slide_ids) {
    slide_row <- cores_to_delete[cores_to_delete$`Slide ID` == slide_id, ]
    
    # Checks if the results from that image should be incorporated or not
    if (nrow(slide_row) == 0) {
      warning(paste("Slide ID", slide_id, "not found in cores_to_delete file. Skipping exclusions for this slide."))
      excluded_cores <- character(0)
    } else {
      excluded_cores <- slide_row %>%
        select(starts_with("cores_to_delete_")) %>%
        unlist() %>%
        unique() %>%
        na.omit() %>%
        as.character()
    }
    
    core_list <- list()
    
    for (core_name in core_names[[tissue]]) {
      if (core_name %in% excluded_cores) {
        next
      }
      
      file_name <- paste0(tissue, "_", core_name, "_", slide_id, ".csv")
      file_path <- file.path(folder_path, file_name)
      
      if (file.exists(file_path)) {
        core <- read.csv(file_path)
        if (!"ROI" %in% colnames(core)) {
          core$ROI <- rep(c("Cytoplasm", "Rootlet", "TZ", "Cilium", "Nucleus"), length.out = nrow(core))
        }
        core_list[[core_name]] <- core
      }
    }
    
    # Takes the median value for each protein and tissue
    if (length(core_list) > 0) {
      combined_core <- bind_rows(core_list)
      final_result_pixel <- combined_core %>%
        group_by(ROI) %>%
        summarize(Median_Category = median(Mean, na.rm = TRUE))
      
      result_list_pixel[[paste(tissue, slide_id, sep = "_")]] <- final_result_pixel
    }
  }
}



# Organize results in list
# Initialize a list to store the results per tissue
Results_per_tissue_PIXEL <- list()

# Loop through each tissue
for (tissue in tissues) {
  # Create a list to store results for this tissue
  tissue_results <- list()
  
  # Loop through each gene (sorterade alfabetiskt)
  sorted_genes <- sort(SlideIDs)
  
  # Loop through each sorted gene
  for (gene in sorted_genes) {
    # Combine the tissue and gene name to access the result
    key <- paste(tissue, gene, sep = "_")
    
    # Check if the result exists for this combination
    if (key %in% names(result_list_pixel)) {
      # Add the result to the list for this tissue
      tissue_results[[gene]] <- result_list_pixel[[key]]
    }
  }
  
  # Add the list of results for this tissue to the nested result list
  Results_per_tissue_PIXEL[[tissue]] <- tissue_results
}


# Prep data for heatmap and clustering
# Extract for each tissue
Fallopian_PIXEL <- Results_per_tissue_PIXEL$Fallopian
Cervix_PIXEL <- Results_per_tissue_PIXEL$Cervix
Endometrium_PIXEL <- Results_per_tissue_PIXEL$Endometrium
Nasopharynx_PIXEL <- Results_per_tissue_PIXEL$Nasopharynx
Bronchus_PIXEL <- Results_per_tissue_PIXEL$Bronchus


# Convert the list of tibbles into a single data frame
combined_df_FT <- bind_rows(
  lapply(names(Fallopian_PIXEL), function(protein) {
    Fallopian_PIXEL[[protein]] %>%
      mutate(Protein = protein)
  }),
  .id = NULL
)
combined_df_C <- bind_rows(
  lapply(names(Cervix_PIXEL), function(protein) {
    Cervix_PIXEL[[protein]] %>%
      mutate(Protein = protein)
  }),
  .id = NULL
)
combined_df_E <- bind_rows(
  lapply(names(Endometrium_PIXEL), function(protein) {
    Endometrium_PIXEL[[protein]] %>%
      mutate(Protein = protein)
  }),
  .id = NULL
)
combined_df_N <- bind_rows(
  lapply(names(Nasopharynx_PIXEL), function(protein) {
    Nasopharynx_PIXEL[[protein]] %>%
      mutate(Protein = protein)
  }),
  .id = NULL
)
combined_df_B <- bind_rows(
  lapply(names(Bronchus_PIXEL), function(protein) {
    Bronchus_PIXEL[[protein]] %>%
      mutate(Protein = protein)
  }),
  .id = NULL
)


# Pivot to a wide format 
matrix_df_FT <- combined_df_FT %>%
  pivot_wider(names_from = ROI, values_from = Median_Category, values_fill = 0) %>%
  rename("Slide ID" = "Protein") # 204
matrix_df_C <- combined_df_C %>%
  pivot_wider(names_from = ROI, values_from = Median_Category, values_fill = 0) %>%
  rename("Slide ID" = "Protein") # 200
matrix_df_E <- combined_df_E %>%
  pivot_wider(names_from = ROI, values_from = Median_Category, values_fill = 0) %>%
  rename("Slide ID" = "Protein") # 192
matrix_df_N <- combined_df_N %>%
  pivot_wider(names_from = ROI, values_from = Median_Category, values_fill = 0) %>%
  rename("Slide ID" = "Protein") # 207
matrix_df_B <- combined_df_B %>%
  pivot_wider(names_from = ROI, values_from = Median_Category, values_fill = 0) %>%
  rename("Slide ID" = "Protein") # 206



# Read in log file from image analysis
# We do this to know which images failed thresholding

# Combine with meta data files
# read in log file from image analysis
log_data <- read_csv("file_here", col_names = FALSE) %>%
  rename("Log" = "X1") 


# Sort out only the rows we need and fix into nice format
# Reference string
reference_string <- "Skipping image: 20466144_Scan1.er_Core[1,8,A]_[10730,55051]_component_data" # This tells us which images it skipped tue to failed thresholding
reference_length <- nchar(reference_string) # Get its length


# This long pipe fixes the whole log to a readable format :)
log_data_clean <- log_data %>% # failed roi/core image by image
  mutate(Log = as.character(Log)) %>% filter(!str_detect(Log, "numLevels|threshold|Warning|Processing|7,A|7,B|8,A|8,B")) %>%
  filter(nchar(Log) <= reference_length) %>%
  mutate(
    is_target = str_detect(Log, "No ROI created for")
  ) %>%
  mutate(
    keep_row = is_target | lead(is_target, 1) | lead(is_target, 2)
  ) %>%
  filter(keep_row) %>%  # Keep only the marked rows
  select(-is_target, -keep_row) %>% # Clean up helper columns
  mutate(
    Log = if_else(
      str_detect(Log, "^[0-9]+_"),                  # Check if row starts with a number followed by an underscore
      str_extract(Log, "^[0-9]+"),                 # Extract only the leading number
      Log                                         # Keep the original value if it doesn't match
    )
  ) %>% 
  mutate(
    Log = if_else(
      str_detect(Log, ".*_.*"),                            # Check if the row contains an underscore
      str_extract(Log, "[^_]+$"),                         # Extract everything after the last underscore
      Log                                                # Keep original value if no underscore is present
    )
  ) %>%
  mutate(
    Log = if_else(
      str_detect(Log, "No ROI created for .*\\."),       # Match rows with this pattern
      str_extract(Log, "(?<=No ROI created for )[a-zA-Z]+(?=\\.)"), # Extract the word after "No ROI created for " and before the dot
      Log                                                # Keep the original value for other rows
    )
  ) %>%
  mutate(New_Column = rep(c("Protein", "Core", "ROI"), length.out = n())) %>%  # Repeat the sequence
  relocate(New_Column, .before = 1) %>%
  mutate(group = cumsum(New_Column == "Protein")) %>%  # Create a group for each Protein-Core-ROI triplet
  pivot_wider(names_from = New_Column, values_from = Log) %>%  # Reshape to wide format
  select(-group)  # Remove the temporary group column


# Now we summarize it a bit more
# make a summary table to have like cores to delete table
# summarize all cores per protein, like 1A,1B,1C...

# Summarize the data
failed_image_analysis <- log_data_clean %>%
  group_by(Protein) %>%  # Group by Protein
  summarize(Cores = paste(Core, collapse = ","), .groups = "drop") %>% rename("Slide ID" = "Protein")


# Combine with cores to delete table (the table with manually removed images during QC just to get a view of how mich data we actually have for each protein)
All_failed_images <- failed_image_analysis %>% right_join(MetaData.Cilia, by = "Slide ID") %>%
  filter(!State %in% c("Failed", "Suggested")) # 216




# Combine all pixel value matrices with this to get the gene name and not only slide ID and to get all information in one place
matrix_FT <- matrix_df_FT %>% inner_join(All_failed_images, by = "Slide ID") %>%
  rename("failed_images" = "Cores")

matrix_E <- matrix_df_E %>% inner_join(All_failed_images, by = "Slide ID") %>%
  rename("failed_images" = "Cores")

matrix_C <- matrix_df_C %>% inner_join(All_failed_images, by = "Slide ID") %>%
  rename("failed_images" = "Cores")

matrix_N <- matrix_df_N %>% inner_join(All_failed_images, by = "Slide ID") %>%
  rename("failed_images" = "Cores")

matrix_B <- matrix_df_B %>% inner_join(All_failed_images, by = "Slide ID") %>%
  rename("failed_images" = "Cores")


# Matrix all tissues 
# Now when we have all information for all proteins in all tissues we want to combine it into one big matrix to generate the heatmap later on
# Make each matrix nice and clean, and rename columns to track which value belongs to which tissue
mx_FT <- matrix_FT %>% select(Gene, Cilium, TZ, Rootlet, Cytoplasm, Nucleus) %>%
  rename("Cilium FT" = "Cilium",
         "TZ FT" = "TZ",
         "Rootlet FT" = "Rootlet",
         "Cytoplasm FT" = "Cytoplasm",
         "Nucleus FT" = "Nucleus")
mx_E <- matrix_E %>% select(Gene, Cilium, TZ, Rootlet, Cytoplasm, Nucleus) %>%
  rename("Cilium E" = "Cilium",
         "TZ E" = "TZ",
         "Rootlet E" = "Rootlet",
         "Cytoplasm E" = "Cytoplasm",
         "Nucleus E" = "Nucleus")
mx_C <- matrix_C %>% select(Gene, Cilium, TZ, Rootlet, Cytoplasm, Nucleus) %>%
  rename("Cilium C" = "Cilium",
         "TZ C" = "TZ",
         "Rootlet C" = "Rootlet",
         "Cytoplasm C" = "Cytoplasm",
         "Nucleus C" = "Nucleus")
mx_N <- matrix_N %>% select(Gene, Cilium, TZ, Rootlet, Cytoplasm, Nucleus) %>%
  rename("Cilium N" = "Cilium",
         "TZ N" = "TZ",
         "Rootlet N" = "Rootlet",
         "Cytoplasm N" = "Cytoplasm",
         "Nucleus N" = "Nucleus")
mx_B <- matrix_B %>% select(Gene, Cilium, TZ, Rootlet, Cytoplasm, Nucleus) %>%
  rename("Cilium B" = "Cilium",
         "TZ B" = "TZ",
         "Rootlet B" = "Rootlet",
         "Cytoplasm B" = "Cytoplasm",
         "Nucleus B" = "Nucleus")



# List of all tables
Matrix_tables <- list(mx_FT, mx_C, mx_E, mx_N, mx_B)

# Perform a full join to get all values into the same matrix
combined_matrix_table <- Reduce(function(x, y) full_join(x, y, by = "Gene"), Matrix_tables) %>% # 190
  drop_na() # 187

# Now we have one big matrix with all image analysis results (median) for each protein
# To be able to put it all into one figure, we need to scale the data


### Scaling
# Load the Excel file
df <- combined_matrix_table


# Separate gene names
genes <- df[, 1]  # First column contains gene names
data_values <- df[, -1]  # All other columns contain numerical values

# Min-Max normalization function applied per row
row_min_max_norm <- function(row) {
  return ((row - min(row, na.rm = TRUE)) / (max(row, na.rm = TRUE) - min(row, na.rm = TRUE)))
}


# We normalize by row to do it for the whole protein to have the difference in intensity
# All tissues are on the same slide stained with the same dilution for the candidate protein, therefore we normalize by row, so by slide one could say
# Apply row-wise normalization
scaled_data <- t(apply(data_values, 1, row_min_max_norm))  # Transpose after applying row-wise function

# Convert back to data frame and keep column names
scaled_data <- as.data.frame(scaled_data)
colnames(scaled_data) <- colnames(data_values)

# Combine gene names with scaled data
scaled_df <- cbind(genes, scaled_data)
# Now we have our big matrix ready for visualization :) 


### Heatmap time!
# Reorder the columns in the heatmap data frame
heatmap_data_ALL <- scaled_df %>% 
  column_to_rownames("Gene") 

# Plot a heatmap quick and easy to get an idea of how it will look
Heatmap.ALL <- pheatmap(heatmap_data_ALL,
                       cluster_cols = T, show_rownames = F)

# Applu clustering by row (stored in new column)
# We chose 7 clusters
ALL.clust <- cbind(heatmap_data_ALL,
                   cluster = cutree(Heatmap.ALL$tree_row,
                                    k = 7))
# Make it into data frame
ALL.clust.df <- as.data.frame(cbind(heatmap_data_ALL,
                                    cluster = cutree(Heatmap.ALL$tree_row,
                                                     k = 7)))

# Re-order original data
# Save row-order
order_ALL <- Heatmap.ALL$tree_row[["order"]]

# Reorder the row names of MX_spt_matrix
ordered_row_names_ALL <- rownames(heatmap_data_ALL)[order_ALL]

# Reorder the rows of ALL.clust.df based on the ordered row names
ALL.clust.df <- ALL.clust.df[match(ordered_row_names_ALL, rownames(ALL.clust.df)), ]

# Plot cluster dendogram
sort(cutree(Heatmap.ALL$tree_row, k=7))
plot(Heatmap.ALL$tree_row)

# Make clusters to factors
ALL.clust.df$cluster <- as.factor(ALL.clust.df$cluster)
ALL.clust.cluster <- ALL.clust.df %>% select(cluster) %>% rename("Cluster" = "cluster") # We rename cluster here to be able to choose colors of cluster bar!

# Define custom colors for HM and clusters
my_colors <- colorRampPalette(c("#ffffff","#d1b3c4", "#b392ac", "#735d78"))(100) # heatmap colors
ann_colors <- list(Cluster = c("1" = "#c0dfa1", "2" = "#f08080", "3" = "#fcbf6e", "4" = "#90c2e7", "5" = "#c19ee0", "6" = "#fff3b0", "7" = "#f4acb7")) # cluster colors

# Ensure column name is correct (avoid duplicate annotations) 
colnames(ALL.clust.cluster) <- "Cluster"

# Generate heatmap 
Heatmap.ALL.col <- pheatmap(
  heatmap_data_ALL, 
  cluster_cols = TRUE, 
  show_rownames = F,  
  border_color = NA,
  annotation_row = ALL.clust.cluster,  # Our clustering
  annotation_colors = ann_colors,  
  main = "All 5 Tissues Scaled row-wise (per slide)", 
  color = my_colors
)

# Now we have a heatmap :) 

# Save clusters for Gene ontology enrichment
# The only thing left to do now is save the clusters to use for the Gene ontology enrichment
ALL.cluster.GO <- ALL.clust.cluster %>% rownames_to_column("Protein") %>%
  mutate(cluster = recode(Cluster, "1" = "Cluster.1",
                          "2" = "Cluster.2",
                          "3" = "Cluster.3",
                          "4" = "Cluster.4",
                          "5" = "Cluster.5",
                          "6" = "Cluster.6",
                          "7" = "Cluster.7")) 

saveRDS("ALL Cluster GO 7 CLUSTERS 2025-02-13 MEDIAN.rds")

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


# End of script :)


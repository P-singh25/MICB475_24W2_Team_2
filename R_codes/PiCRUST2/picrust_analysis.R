##### Install Necessary Packages Once #####

# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Create a list of all the packages needed to install and then install them 
pkgs <- c("ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

for (pkg in pkgs) {
   if (!requireNamespace(pkg, quietly = TRUE))
       BiocManager::install(pkg)
}
# when asked if you want to update all, some or none please type "n" for none

# After installing all of its above dependencies, Install ggpicrust2
install.packages("ggpicrust2")



##### Load packages #####


# Load all necessary libraries
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)
library(dplyr)
library(ggplot2)


##### Import files and preparing tables #####


#Importing the pathway PICrsut2
abundance_file <- "PiCRUST2/pathway_abundance.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE, skip = 1)
abundance_data  =as.data.frame(abundance_data)

#Import the metadata file
metadata <- read_delim("PiCRUST2/melanoma_metadata_updated_copy.tsv")

#Change the name of the column to sample-id in meta data and #OTU-id in abundance data
colnames(abundance_data)[colnames(abundance_data) == "#OTU ID"] <- "pathway"
colnames(metadata)[colnames(metadata) == "sample-id"] <- "sample_id"

# Filter the metadata for only 'Pre-ICI' and 'Post-ICI3' groups
filtered_metadata <- metadata %>%
  filter(group %in% c("Pre-ICI", "Post-ICI3")) %>%
  filter(location == "Spleen")

# Save the filtered metadata to as a new file
write_delim(filtered_metadata, "PiCRUST2/filtered_metadata_for_picrust.tsv", delim = "\t")

#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names = filtered_metadata$'sample_id'
sample_names = append(sample_names, "pathway")
abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered =  abundance_data_filtered[, colSums(abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered) = NULL

#Verify samples in metadata match samples in abundance_data
abun_samples = rownames(t(abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
metadata_final = filtered_metadata[filtered_metadata$`sample_id` %in% abun_samples,] #making sure the filtered metadata only includes these samples




#### DESEq ####


#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), 
                                        metadata = metadata_final, group = "group", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                       daa_results_df = abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
feature_with_p_0.05 <- abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
feature_desc = inner_join(feature_with_p_0.05,metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% 
  filter(`pathway` %in% feature_with_p_0.05$feature)  
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description
#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
abundance_desc = abundance_desc[,-c(85:ncol(abundance_desc))] 



#### Heatmap ####

# Get just the sample columns + 'feature'
sample_cols <- c("feature", metadata_final$sample_id)
abundance_desc_filtered <- abundance_desc_[, colnames(abundance_desc) %in% sample_cols]

# Generate a heatmap and save it
pathway_heatmap(abundance = abundance_desc_filtered %>% column_to_rownames("feature"), metadata = metadata_final, group = "group")

ggsave("../pathway_heatmap.png", 
       plot = last_plot() + 
         theme(axis.text.y = element_text(size = 5, hjust = 1)),  
       width = 12, height = 10, dpi = 300)

# Reduce Number of Pathways, keep only significant
top_pathways = abundance_desc_filtered %>%
  slice_max(order_by = rowSums(abs(.[,-1])), n = 20) # Select top 20

# Re-run pathway_heatmap() with filtered data and save it 
heatmap_plot <- pathway_heatmap(abundance = top_pathways %>% column_to_rownames("feature"), 
                metadata = metadata_final, group = "group")

# Modify the plot to use red gradient only (from white to red)
heatmap_red_only <- heatmap_plot + 
  scale_fill_gradient(
    low  = "white",
    high = "red"
  )

ggsave("../pathway_heatmap_30.png", 
       plot = last_plot() + 
         theme(axis.text.y = element_text(size = 9, hjust = 1)),  
       width = 12, height = 10, dpi = 300)




#### PCA plot ####


# Create a new data frame of filtered abundance data with no name column name 
abundance_pca <- abundance_data_filtered %>% column_to_rownames("pathway")

# Create a new data frame of metadata where the column name 'sample_id' is changed to 'sample_name' 
metadata_pca <- metadata_final
colnames(metadata_pca)[colnames(metadata_pca) == "sample_id"] <- "sample_name"

# Remove rows (pathways) with zero variance across all samples
abundance_pca_clean <- abundance_pca[apply(abundance_pca, 1, function(x) var(x) != 0), ]

# Generate pathway PCA plot and save it 
pathway_pca(abundance = abundance_pca_clean, metadata = metadata_pca, group = "group")

ggsave("../pathway_pca_plot.png", width = 8, height = 6, dpi = 300)



#### Bar plot ####


# Generating a bar plot representing log2FC from the custom deseq2 function
# In the Deseq2 function script and the metadata category of interest which is group has been updated 

# Lead the function in
source("Picrust2/DESeq2_function.R")

# Run the function on your own data
res =  DEseq2_function(abundance_data_filtered, metadata_final, "group")
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]
View(res_desc)

# Filter to only include significant pathways
sig_res = res_desc %>%
  filter(pvalue < 0.05)
  
# You can also filter by Log2fold change
# keep only greater than 2 fold change or less than -2
res_desc$pvalue <- as.numeric(res_desc$pvalue)
res_desc$log2FoldChange <- as.numeric(res_desc$log2FoldChange)
sig_res = res_desc %>%
  filter(pvalue < 0.05 & (log2FoldChange > 2 | log2FoldChange < -2 ))

# Generate the bar plot for our data and save it 
sig_res <- sig_res[order(sig_res$log2FoldChange),]
bar <- ggplot(data = sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

ggsave("../pathway_barplot.png", plot = bar, width = 8, height = 6, dpi = 300)
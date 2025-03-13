# Load required packages
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tibble)

#### Taxonomy bar plot ####

# Load the processed phyloseq object
load("pj2.RData")

# Remove Unclassified Taxa & Aggregate Data at Phylum Level
pj2_rare <- pj2 %>%
  subset_taxa(!apply(is.na(tax_table(pj2)) | 
                       tax_table(pj2) == "" | 
                       tax_table(pj2) == "Unclassified" | 
                       tax_table(pj2) == "Unassigned", 1, all)) %>%
  tax_glom(taxrank = "Phylum", NArm = FALSE)  # Collapse to Phylum level

# Create a Combined "group_location" Variable Before Merging
sample_data(pj2_rare)$group_location <- paste(
  sample_data(pj2_rare)$group, 
  sample_data(pj2_rare)$location, 
  sep = "_"
)

# Merge Samples by "group_location"
pj2_grouped <- merge_samples(pj2_rare, "group_location")

# Restore 'group' and 'location' Separately
sample_data(pj2_grouped)$group <- gsub("_.*", "", rownames(sample_data(pj2_grouped)))  # Extract group
sample_data(pj2_grouped)$location <- gsub(".*_", "", rownames(sample_data(pj2_grouped)))  # Extract location

## Relative Abundance ## 


# Convert to Relative Abundance
pj2_RA_grouped <- transform_sample_counts(pj2_grouped, function(x) x / sum(x))

# Select Top 10 Most Abundant Phyla and Collapse Others into "Other"
# Calculate mean relative abundance of each phylum
phyla_means <- taxa_sums(pj2_RA_grouped) / nsamples(pj2_RA_grouped)
top10_phyla <- names(sort(phyla_means, decreasing = TRUE)[1:10])  # Select Top 10 Phyla

# Collapse the remaining phyla into "Other"
tax_table(pj2_RA_grouped)[!rownames(tax_table(pj2_RA_grouped)) %in% top10_phyla, "Phylum"] <- "Other"

# Ensure 'group' is a Factor with Correct Order
sample_data(pj2_RA_grouped)$group <- factor(
  sample_data(pj2_RA_grouped)$group, 
  levels = c("Day 0", "Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3")
)

# Generate the Taxonomic Bar Plot
taxa_relative <- plot_bar(pj2_RA_grouped, fill = "Phylum", x = "group") +  
  theme_minimal() +
  facet_wrap(.~location, scales = "free") +  
  labs(
    title = "Taxonomic Bar Plot (Relative Abundance)", 
    x = "ICI Treatment Group", 
    y = "Relative Abundance", 
    fill = "Phylum"
  ) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = NA), # To have white background in the whole plot
    plot.background = element_rect(fill = "white", color = NA),  
    legend.background = element_rect(fill = "white"),  
    strip.background = element_rect(fill = "white"),  
  ) +  
  scale_x_discrete(labels = c("Day 0", "Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3")) +  
  scale_y_continuous(labels = scales::percent_format())

# Print the relative abundance bar plot
print(taxa_relative )

# Save the plot
ggsave("taxonomy_relative_abundance.png", taxa_relative , height=8, width=12)





#### Taxonomy 1 -- failed ####

# Ensure 'group' shows up in the correct order
sample_data(pj2_rare)$group <- factor(sample_data(pj2_rare)$group, 
                                      levels = c("Day 0", "Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3"))

# Convert to relative abundance
pj2_RA <- transform_sample_counts(pj2_rare, function(x) x / sum(x))


# Plot bar plot of taxonomy (Relative Abundance)
gg_taxa_relative <- plot_bar(pj2_RA, x="group", fill="Phylum") + 
  facet_wrap(.~location, scales = "free") +  
  theme_minimal() +
  labs(title = "Taxonomic Bar Plot (Relative Abundance)",
       x = "ICI Treatment Group",
       y = "Relative Abundance",
       fill = "Phylum") +
  theme(
    panel.background = element_rect(fill = "white", color = NA), # To have white background in the whole plot
    plot.background = element_rect(fill = "white", color = NA),  
    legend.background = element_rect(fill = "white"),  
    strip.background = element_rect(fill = "white"),  
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),  # To remove major and minor grid lines
    panel.grid.minor = element_blank()  
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))  # Ensure bars fill the plot fully

# Print the relative abundance bar plot
print(gg_taxa_relative)

# Save the plot
ggsave("taxonomy_relative_abundance.png", gg_taxa_relative, height=8, width=12)







#### Taxonomy 2 -- failed ####
otu_absolute <- as.data.frame(otu_table(pj2))

phyloseq_tax <- tax_table(pj2)
tax_df <- as.data.frame(phyloseq_tax)
tax_df1 <- rownames_to_column(tax_df, var = "ASV")

phyloseq_sam <- sample_data(pj2)
phyloseq_sam_df <- as.data.frame(phyloseq_sam)

otu_absolute1 <- rownames_to_column(otu_absolute, var = "ASV")

otu_absolute_melt <- reshape2::melt(data = otu_absolute1, 
                                    measure.vars = 2:ncol(otu_absolute1),  # All sample columns
                                    variable.name = "Sample", 
                                    value.name = "Absolute_Abundance", 
                                    as.is = TRUE)

taxa_ab <- dplyr::left_join(otu_absolute_melt, tax_df1, by = "ASV")
phyloseq_sam_sub <- phyloseq_sam_df[, c("location", "group")]  # Replace with relevant column names
phyloseq_sam_sub1 <- rownames_to_column(phyloseq_sam_sub, var = "Sample")
taxa_ab_full <- dplyr::left_join(taxa_ab, phyloseq_sam_sub1, by = "Sample")

taxa_plot <- ggplot(data = taxa_ab_full, aes(x = group, y = Absolute_Abundance)) +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +  # Stacked bar plot
  ylab('Absolute abundance (read counts)') +  # Y-axis label
  facet_wrap(~ location, scales = "free_x") +  # Facet by location
  ggtitle("Taxa Bar Plot by Location") +  # Plot title
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # Rotate x-axis labels for clarity
  theme_bw()

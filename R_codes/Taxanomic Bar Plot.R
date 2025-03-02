# Load required packages
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyverse)

# Load the processed phyloseq object
load("pj2.RData")

# Using the rarefied dataset, filter out completely unclassified taxas and then aggregate the data at the Phylum level
pj2_rare <- rare %>%
  subset_taxa(!apply(is.na(tax_table(rare)) | tax_table(rare) == "" | tax_table(rare) == "Unclassified"| tax_table(rare) == "Unassigned", 1, all)) %>%
  tax_glom(taxrank = "Phylum", NArm = FALSE)

#### Taxonomy bar plots ####

## Absolute Abundance ##

# To ensure 'group' shows up in the correct order
sample_data(pj2_rare)$group <- factor(sample_data(pj2_rare)$group, 
                                      levels = c("Day 0", "Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3"))

# Plot bar plot of taxonomy (Absolute Abundance)
gg_taxa_absolute <- plot_bar(pj2_rare, x="group", fill="Phylum") + 
  facet_wrap(.~location, scales = "free") +  
  theme_minimal() +
  labs(title = "Taxonomic Bar Plot (Absolute Abundance)",
       x = "ICI Treatment Group",
       y = "Absolute Abundance",
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
  
# Print the absolute abundance bar plot
print(gg_taxa_absolute)

# Save the plot
ggsave("taxonomy_absolute_abundance.png", gg_taxa_absolute, height=8, width=12)


## Relative Abundance ## - to compare to the paper

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

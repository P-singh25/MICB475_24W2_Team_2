# Load required packages
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tibble)

# Load the processed phyloseq object
load("pj2.RData")
otu_table <- otu_table(pj2)
head(otu_table)
taxa_are_rows(otu_table)

summary(as.vector(otu_table))
sum(rowSums(otu_table) == 0)

sample_sums(pj2)

sample_data(pj2)
otu_table(pj2)



ps_phylum <- tax_glom(ps_relabund, taxrank = "Phylum")

tax_table(pj2_new)
sample_data(pj2_new)

sum(is.na(tax_table(pj2)[, "Phylum"]))
pj2_new <- subset_taxa(pj2, !is.na(Phylum)) 

sum(is.na(tax_table(pj2_new)[, "Phylum"]))

plot_bar(pj2_new, fill="Phylum") 
pj2_RA <- transform_sample_counts(pj2_2, function(x) x/sum(x))

plot_bar(pj2_RA, x = "group", fill="Phylum") + 
  facet_wrap(.~location, scale = "free_x")




# Using the rarefied dataset, filter out completely unclassified taxas and then aggregate the data at the Phylum level
pj2_rare <- pj2 %>%
  subset_taxa(!apply(is.na(tax_table(pj2)) | tax_table(pj2) == "" | tax_table(pj2) == "Unclassified"| tax_table(pj2) == "Unassigned", 1, all)) %>%
  tax_glom(taxrank = "Phylum", NArm = FALSE)

#### Taxonomy bar plots ####

## Absolute Abundance ##

# To ensure 'group' shows up in the correct order
sample_data(pj2_rare)$group <- factor(sample_data(pj2_rare)$group, 
                                      levels = c("Day 0", "Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3"))

# Plot bar plot of taxonomy (Absolute Abundance)
gg_taxa_absolute <- plot_bar(pj2_rare, x="group", fill="Phylum") + 
  facet_wrap(.~location, scales = "free_x") + 
  labs(title = "Taxonomic Bar Plot (Absolute Abundance)",
       x = "ICI Treatment Group",
       y = "Absolute Abundance",
       fill = "Phylum") +
  theme_minimal() +
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

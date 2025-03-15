# Load necessary libraries
library(phyloseq)
library(tidyverse)

load("pj2.RData")

# Step 1: Remove Unclassified Taxa & Aggregate Data at Phylum Level

pj2_rare <- pj2 %>%
  subset_taxa(!apply(is.na(tax_table(pj2)) | 
                       tax_table(pj2) == "" | 
                       tax_table(pj2) == "Unclassified" | 
                       tax_table(pj2) == "Unassigned", 1, all)) %>%
  tax_glom(taxrank = "Phylum", NArm = FALSE)  # Collapse to Phylum level

pj2_filtered <- subset_samples(pj2_rare, location != "Stool" & group != "Day 0")


# Step 2: Create a Combined "group_location" Variable Before Merging
sample_data(pj2_filtered)$group_location <- paste(
  sample_data(pj2_filtered)$group, 
  sample_data(pj2_filtered)$location, 
  sep = "_"
)

# Step 3: Merge Samples by "group_location"
pj2_grouped <- merge_samples(pj2_filtered, "group_location")

# Step 4: Restore 'group' and 'location' Separately
sample_data(pj2_grouped)$group <- gsub("_.*", "", rownames(sample_data(pj2_grouped)))  # Extract group
sample_data(pj2_grouped)$location <- gsub(".*_", "", rownames(sample_data(pj2_grouped)))  # Extract location

# Step 5: Convert to Relative Abundance
pj2_RA_grouped <- transform_sample_counts(pj2_grouped, function(x) x / sum(x))

# Step 6: Select Top 10 Most Abundant Phyla and Collapse Others into "Other"
# Calculate mean relative abundance of each phylum
phyla_means <- taxa_sums(pj2_RA_grouped) / nsamples(pj2_RA_grouped)
top10_phyla <- names(sort(phyla_means, decreasing = TRUE)[1:10])  # Select Top 10 Phyla

# Collapse the remaining phyla into "Other"
tax_table(pj2_RA_grouped)[!rownames(tax_table(pj2_RA_grouped)) %in% top10_phyla, "Phylum"] <- "Other"

# Step 7: Ensure 'group' is a Factor with Correct Order
sample_data(pj2_RA_grouped)$group <- factor(
  sample_data(pj2_RA_grouped)$group, 
  levels = c("Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3")
)

sample_data(pj2_RA_grouped)$location <- factor(sample_data(pj2_RA_grouped)$location, levels = c("Spleen", "Tumor","TDLN", "MLN"))


# Step 8: Generate the Taxonomic Bar Plot
taxa_bar_plot <- plot_bar(pj2_RA_grouped, fill = "Phylum", x = "group") +  
  theme_classic() +
  facet_wrap(.~location, scales = "free", nrow = 1) +  
  labs(
    # title = "Taxonomic Bar Plot (Relative Abundance)", 
    x = "ICI Treatment Group", 
    y = "Relative Abundance", 
    fill = "Phylum"
  ) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.title = element_text(size = 14, hjust = 0.5),
    legend.text = element_text(size = 12), 
    axis.line = element_line(size = 0),
    strip.background = element_rect(color = "white", fill = "white", size = 1), 
    panel.border = element_rect(color = "black", fill = NA)
  ) +  
  # scale_x_discrete(labels = c("Day 0", "Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3")) +  
  scale_y_continuous(labels = scales::percent_format())

### reorder the phylum

# Convert phyloseq object to a dataframe
df <- psmelt(pj2_RA_grouped)

# Reorder Phylum levels
df$Phylum <- factor(df$Phylum, levels = c(setdiff(unique(df$Phylum), "Other"), "Other"))

# Plot using the updated dataframe
taxa_bar_plot <- ggplot(df, aes(x = group, y = Abundance, fill = Phylum)) +  
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  facet_wrap(.~location, scales = "free", nrow = 1) +  
  labs(
    x = "ICI Treatment Group", 
    y = "Relative Abundance", 
    fill = "Phylum"
  ) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.title = element_text(size = 14, hjust = 0.5),
    legend.text = element_text(size = 12), 
    axis.line = element_line(size = 0),
    strip.background = element_rect(color = "white", fill = "white", size = 1), 
    panel.border = element_rect(color = "black", fill = NA)
  ) +  
  scale_y_continuous(labels = scales::percent_format())

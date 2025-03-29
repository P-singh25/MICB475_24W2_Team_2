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
df <- psmelt(pj2_RA_grouped)

# Reorder Phylum levels
df$Phylum <- factor(df$Phylum, levels = c(setdiff(unique(df$Phylum), "Other"), "Other"))

taxa_bar_plot_1 <- plot_bar(pj2_RA_grouped, fill = "Phylum", x = "group") +  
  theme_bw() +
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
    # strip.background = element_rect(color = "white", fill = "white", size = 1), 
    # panel.border = element_rect(color = "black", fill = NA)
  ) +  
 scale_fill_manual(
  values = c(
   "p__Bacteroidota" = "#1b7838",
  "p__Firmicutes" = "#a6dca1",
 "p__Proteobacteria" = "#fe9f9a",
 "p__Actinobacteriota" = "#f3694e",
 "p__Verrucomicrobiota" = "#dd65af",
 "p__Fusobacteriota" = "#9281ff",
 "p__Deinococcota" = "#0e9ee2",
 "p__Euryarchaeota" = "#00d69f",
 "p__Planctomycetota" = "#ffd165",
 "p__Cyanobacteria" = "#1b9cd6",
 "Other" = "gray70"
 )
  ) +
  # scale_x_discrete(labels = c("Day 0", "Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3")) +
  scale_y_continuous(labels = scales::percent_format())




df <- psmelt(pj2_RA_grouped)
unique(df$Phylum)






### reorder the phylum

# Convert phyloseq object to a dataframe
df <- psmelt(pj2_RA_grouped)

# Reorder Phylum levels
df$Phylum <- factor(df$Phylum, levels = c(setdiff(unique(df$Phylum), "Other"), "Other"))
df$Phylum <- gsub("p__", "", df$Phylum)
df$Phylum <- factor(df$Phylum, levels = c(
   "Firmicutes", "Proteobacteria", "Bacteroidota", "Actinobacteriota",
  "Verrucomicrobiota", "Fusobacteriota", "Deinococcota", "Euryarchaeota",
  "Planctomycetota", "Cyanobacteria", "Other"
))

# Plot using the updated dataframe
taxa_bar_plot <- ggplot(df, aes(x = group, y = Abundance, fill = Phylum)) +  
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  facet_wrap(.~location, scales = "free", nrow = 1) +  
  labs(
    x = "ICI Treatment Group", 
    y = "Relative Abundance", 
    fill = "Phylum"
  ) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, hjust = 0.5, color = "black"),
    legend.text = element_text(size = 12, color = "black"), 
    # axis.line = element_line(size = 0),
    #strip.background = element_rect(color = "white", fill = "white", size = 1), 
    #panel.border = element_rect(color = "black", fill = NA)
  ) + 
  scale_fill_manual(
    values = c(
     "Bacteroidota" = "#ffd165",
    "Firmicutes" = "#a6dca1",
   "Proteobacteria" = "#fe9f9a",
  "Actinobacteriota" = "#f3694e",
  "Verrucomicrobiota" = "#0e9ee0",
   "Fusobacteriota" = "#00d69f",
  "Deinococcota" = "#9381ff",
   "Euryarchaeota" = "#dd65af",
  "Planctomycetota" = "#fad3e2",
   "Cyanobacteria" = "#2b7db8",
      "Other" = "gray70"
    )
  ) +
  # scale_x_discrete(labels = c("Day 0", "Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3")) +  
  scale_y_continuous(labels = scales::percent_format())  

taxa_bar_plot
ggsave("Taxa_Plot_Final.png")

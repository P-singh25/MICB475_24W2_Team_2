# load library
library(phyloseq)
library(picante)
library(tidyverse)
library(ape)
library(ggsci)
library(scales)

# load the pj2 R data 
load("pj2.RData")
sample_data(pj2)

#### Beta diversity #####

### filter data ###
# remove the stool and D0 sample 
pj2_filtered <- subset_samples(pj2, location != "Stool" & group != "Day 0")

# Reorganize the order of the location 
sample_data(pj2_filtered)$location <- factor(sample_data(pj2_filtered)$location, 
                                          levels = c("Spleen", "Tumor", "TDLN", "MLN"))  
  

####### Bray curtis #######
### filter data ###
pj2_filtered <- subset_samples(pj2, location != "Stool" & group != "Day 0")

sample_data(pj2_filtered)$location <- factor(sample_data(pj2_filtered)$location, 
                                             levels = c("Spleen", "Tumor", "TDLN", "MLN"))  

zero_samples <- sample_sums(pj2_filtered) == 0
sum(zero_samples)  # Count of empty samples
sample_names(pj2_filtered)[zero_samples]  # List of empty samples

pj2_filtered <- prune_samples(sample_sums(pj2_filtered) > 0, pj2_filtered)

bc_dm <- distance(pj2_filtered, method="bray")  # Bray-Curtis distance
pcoa_bray <- ordinate(pj2_filtered, method="PCoA", distance=bc_dm)
plot_ordination(pj2_filtered, pcoa_bray, color = "location", shape="group")

sample_data(pj2_filtered)$group <- factor(sample_data(pj2_filtered)$group, 
                                          levels = c("Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3"))

bcurtis_facet <- plot_ordination(pj2_filtered, pcoa_bray, color = "location") +
  facet_wrap(~ group) +  # Single row layout
  labs(pch = "Treatment group", col = "Location") +
  theme_classic() +
  stat_ellipse(aes(group = location), level = 0.8, linetype = "solid") +
 # coord_fixed(ratio = 1.2) +
  scale_color_npg() +
  theme(legend.title = element_text(size = 14, hjust = 0.5), 
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.y = element_text(size = 12, angle = 0, color = "black"),
        axis.text.x = element_text(size = 12, angle = 0, color = "black"),
        axis.line = element_line(size = 0),
        strip.background = element_rect(color = "white", fill = "white", size = 1),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 12))
bcurtis_facet
ggsave("bray_curtis.png")

# Significance for pre-ici
pj2_group_1 <- subset_samples(pj2_filtered, group == "Pre-ICI")
meta_data_group <- as(sample_data(pj2_group_1),"data.frame")
uu_dm_group <- distance(pj2_group_1, method = "bray")
perm_results_pre <- adonis2(uu_dm_group ~ location, data = meta_data_group, permutations = 999)
print(perm_results_pre)

# Significance for post-ici3
pj2_group <- subset_samples(pj2_filtered, group == "Post-ICI3")
meta_data_group <- as(sample_data(pj2_group),"data.frame")
uu_dm_group <- distance(pj2_group, method = "bray")
perm_results_post3 <- adonis2(uu_dm_group ~ location, data = meta_data_group, permutations = 999)
print(perm_results_post3)

###significance of spleen vs others
# Create a new column in metadata: "Spleen" vs. "Others"
sample_data(pj2_filtered)$location_binary <- ifelse(sample_data(pj2_filtered)$location == "Spleen", "Spleen", "Other")

# Ensure it is a factor
sample_data(pj2_filtered)$location_binary <- factor(sample_data(pj2_filtered)$location_binary, levels = c("Spleen", "Other"))

# Subset Pre-ICI group
pj2_group_1 <- subset_samples(pj2_filtered, group == "Pre-ICI")
meta_data_group <- as(sample_data(pj2_group_1), "data.frame")
uu_dm_group <- distance(pj2_group_1, method = "bray")

# PERMANOVA comparing "Spleen" vs. "Other" locations
perm_results_pre <- adonis2(uu_dm_group ~ location_binary, data = meta_data_group, permutations = 999)
print(perm_results_pre)

# Subset Post-ICI3 group
pj2_group <- subset_samples(pj2_filtered, group == "Post-ICI3")
meta_data_group <- as(sample_data(pj2_group), "data.frame")
uu_dm_group <- distance(pj2_group, method = "bray")

# PERMANOVA comparing "Spleen" vs. "Other" locations
perm_results_post3 <- adonis2(uu_dm_group ~ location_binary, data = meta_data_group, permutations = 999)
print(perm_results_post3)

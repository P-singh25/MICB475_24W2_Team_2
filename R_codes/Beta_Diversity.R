library(phyloseq)
library(picante)
library(tidyverse)
library(ape)
library(ggsci)
library(scales)


load("pj2.RData")
sample_data(pj2)

#### Beta diversity #####

### filter data ###
pj2_filtered <- subset_samples(pj2, location != "Stool" & group != "Day 0")

sample_data(pj2_filtered)$location <- factor(sample_data(pj2_filtered)$location, 
                                          levels = c("Spleen", "Tumor", "TDLN", "MLN"))  
  

#### Unweighted UniFrac ####
# Compute unweighted UniFrac distance
uu_dm <- distance(pj2_filtered, method="unifrac")

# Check for NA or infinite values
sum(is.na(uu_dm))  # Should return 0
sum(is.infinite(uu_dm))  # Should return 0

# Replace NA values if needed
uu_dm[is.na(uu_dm)] <- 0

# Prune samples with zero counts
pj2_filtered <- prune_samples(sample_sums(pj2_filtered) > 0, pj2_filtered)

# Recompute distance matrix after pruning
uu_dm <- distance(pj2_filtered, method="unifrac")

# Perform PCoA ordination with Cailliez correction
pcoa_uu <- ordinate(pj2_filtered, method="PCoA", distance=uu_dm, correction="cailliez")

# Plot ordination
plot_ordination(pj2_filtered, pcoa_uu, color = "location", shape="group")

#make Pre-ICI the first in the plot
sample_data(pj2_filtered)$group <- factor(sample_data(pj2_filtered)$group, 
                                          levels = c("Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3"))

# Faceted plot
# unifrac_facet <- plot_ordination(pj2_filtered, pcoa_uu, color = "location") +
#   facet_wrap(~ group) +
#   labs(pch="Treatment group", col="Organ") +
#   theme_bw() +
#   stat_ellipse(aes(group = location), level = 0.95, linetype = "solid")

unifrac_facet <- plot_ordination(pj2_filtered, pcoa_uu, color = "location") +
  facet_wrap(~ group, nrow = 2, scales = "fixed") +  # Single row layout
  labs(pch = "Treatment group", col = "Location") +
  theme_classic() +
  stat_ellipse(aes(group = location), level = 0.8, linetype = "solid") +
  coord_fixed(ratio = 1.2) +
  scale_color_npg() +
  theme(legend.title = element_text(size = 14, hjust = 0.5), 
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.y = element_text(size = 12, angle = 0, color = "black"),
        axis.text.x = element_text(size = 12, angle = 0, color = "black"),
        axis.line = element_line(size = 0),
        strip.background = element_rect(color = "white", fill = "white", size = 1),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 12))


  

unifrac_facet
ggsave("Beta_Diversity_unweighted_unifrac.png")  # Adjust size if needed


unifrac_facet
ggsave("unifrac_facet.png", width = 10, height = 8)  # Adjust size if needed

unifrac_facet
ggsave("unifrac_facet.png")

# Base plot
unifrac_base <- plot_ordination(pj2_filtered, pcoa_uu, color = "location", shape="group") +
  labs(pch="Treatment group", col="Organ") +
  theme_bw() +
  stat_ellipse(aes(group = location), level = 0.95, linetype = "solid")

unifrac_base

# Significance for post-ici3
pj2_group <- subset_samples(pj2_filtered, group == "Post-ICI3")
meta_data_group <- as(sample_data(pj2_group),"data.frame")
uu_dm_group <- distance(pj2_group, method = "unifrac")
perm_results_post3 <- adonis2(uu_dm_group ~ location, data = meta_data_group, permutations = 999)
print(perm_results_post3)

pj2_group_1 <- subset_samples(pj2_filtered, group == "Pre-ICI")
meta_data_group <- as(sample_data(pj2_group_1),"data.frame")
uu_dm_group <- distance(pj2_group_1, method = "unifrac")
perm_results_pre <- adonis2(uu_dm_group ~ location, data = meta_data_group, permutations = 999)
print(perm_results_pre)

# Weight Unifrac #
wu_dm <- distance(pj2, method="wunifrac")

sum(is.na(wu_dm))  # Should return 0
sum(is.infinite(wu_dm))  # Should return 0

wu_dm[is.na(wu_dm)] <- 0

sum(is.na(wu_dm))

sample_sums(pj2)

pj2 <- prune_samples(sample_sums(pj2) > 0, pj2)

wu_dm <- distance(pj2, method="wunifrac")

pcoa_wu <- ordinate(pj2, method="PCoA", distance=wu_dm, correction="cailliez")

plot_ordination(pj2, pcoa_wu, color = "location", shape="group")

wunifrac_facet <- plot_ordination(pj2, pcoa_wu, color = "location", shape="group") +
  facet_wrap(~ group) +  # Creates a separate plot for each group
  stat_ellipse(aes(group = location), level = 0.95, linetype = "solid") +
  labs(pch="Treatment group", col="Organ") +
  theme_bw() 

wunifrac_facet

wunifrac_base <- plot_ordination(pj2, pcoa_wu, color = "location", shape="group") +
  labs(pch="Treatment group", col="Organ") +
  #stat_ellipse(aes(group = location), level = 0.95, linetype = "solid") +
  theme_bw()

wunifrac_base






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

bcurtis_facet <- plot_ordination(pj2_filtered, pcoa_bray, color = "location", shape = "group") +
  facet_wrap(~ group) +  # Creates separate plots for each group
  #stat_ellipse(aes(group = location), level = 0.95, linetype = "solid") +
  labs(pch="Treatment group", col="Organ") +
  theme_minimal() +
  stat_ellipse(aes(group = location), level = 0.95, linetype = "solid")  # Optional: clean theme

bcurtis_facet

bcurtis_base <- plot_ordination(pj2_filtered, pcoa_bray, color = "location", shape = "group") +
  #stat_ellipse(aes(group = location), level = 0.95, linetype = "solid") +
  labs(pch="Treatment group", col="Organ") +
  theme_minimal()  # Optional: clean theme

bcurtis_base








#unweighted
#### Load Data ####
load("pj2.RData")

### Filter data ###
pj2_filtered <- pj2 %>%
  subset_samples(location != "Stool") %>%
  subset_samples(group != "Day 0")

#### Check Sample Data ####
sample_data(pj2_filtered)
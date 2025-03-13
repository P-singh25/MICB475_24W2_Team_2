library(phyloseq)
library(picante)
library(tidyverse)
library(ape)


load("pj2.RData")
sample_data(pj2)

#### Beta diversity #####

### filter data ###
pj2_filtered <- subset_samples(pj2, location != "Stool" & group != "Day 0")

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

# Faceted plot
unifrac_facet <- plot_ordination(pj2_filtered, pcoa_uu, color = "location", shape="group") +
  facet_wrap(~ group) +
  labs(pch="Treatment group", col="Organ") +
  theme_bw() +
  stat_ellipse(aes(group = location), level = 0.95, linetype = "solid")

unifrac_facet
ggsave("unifrac_facet.png")

# Base plot
unifrac_base <- plot_ordination(pj2_filtered, pcoa_uu, color = "location", shape="group") +
  labs(pch="Treatment group", col="Organ") +
  theme_bw() +
  stat_ellipse(aes(group = location), level = 0.95, linetype = "solid")

unifrac_base







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


# Bray curtis #
zero_samples <- sample_sums(pj2) == 0
sum(zero_samples)  # Count of empty samples
sample_names(pj2)[zero_samples]  # List of empty samples

pj2 <- prune_samples(sample_sums(pj2) > 0, pj2)

bc_dm <- distance(pj2, method="bray")  # Bray-Curtis distance
pcoa_bray <- ordinate(pj2, method="PCoA", distance=bc_dm)
plot_ordination(pj2, pcoa_bray, color = "location", shape="group")

bcurtis_facet <- plot_ordination(pj2, pcoa_bray, color = "location", shape = "group") +
  facet_wrap(~ group) +  # Creates separate plots for each group
  #stat_ellipse(aes(group = location), level = 0.95, linetype = "solid") +
  labs(pch="Treatment group", col="Organ") +
  theme_minimal() +
  stat_ellipse(aes(group = location), level = 0.95, linetype = "solid")  # Optional: clean theme

bcurtis_facet

bcurtis_base <- plot_ordination(pj2, pcoa_bray, color = "location", shape = "group") +
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
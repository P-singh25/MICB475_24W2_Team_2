library(phyloseq)
library(picante)
library(tidyverse)
library(ape)


load("pj2.RData")
sample_data(pj2)

#### Beta diversity #####
# Weight Unifrac #
wu_dm <- distance(pj2, method="wunifrac")

sum(is.na(wu_dm))  # Should return 0
sum(is.infinite(wu_dm))  # Should return 0

wu_dm[is.na(wu_dm)] <- 0

sum(is.na(wu_dm))

sample_sums(pj2)

pj2 <- prune_samples(sample_sums(pj2) > 0, pj2)

pcoa_wu <- ordinate(pj2, method="PCoA", distance=wu_dm, correction="cailliez")

wu_dm <- distance(pj2, method="wunifrac")

plot_ordination(pj2, pcoa_wu, color = "location", shape="group")

gg_pcoa <- plot_ordination(pj2, pcoa_wu, color = "location", shape="group") +
  facet_wrap(~ group) +  # Creates a separate plot for each group
  labs(pch="Treatment group", col="Organ") +
  theme_bw()

gg_pcoa



# Bray curtis #
zero_samples <- sample_sums(pj2) == 0
sum(zero_samples)  # Count of empty samples
sample_names(pj2)[zero_samples]  # List of empty samples

pj2 <- prune_samples(sample_sums(pj2) > 0, pj2)

bc_dm <- distance(pj2, method="bray")  # Bray-Curtis distance
pcoa_bray <- ordinate(pj2, method="PCoA", distance=bc_dm)
plot_ordination(pj2, pcoa_bray, color = "location", shape="group")

gg_pcoa <- plot_ordination(pj2, pcoa_bray, color = "location", shape = "group") +
  facet_wrap(~ group) +  # Creates separate plots for each group
  theme_minimal()  # Optional: clean theme

gg_pcoa



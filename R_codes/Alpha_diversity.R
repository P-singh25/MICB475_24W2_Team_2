library(phyloseq)
library(picante)
library(tidyverse)

#### load phyloseq object ####
load("pj2.RData")
ls()
sample_data(pj2)

#### shannon diversity ####
alpha_diversity <- estimate_richness(pj2, measures = c("Shannon", "Simpson", "Chao1", "Observed"))
alpha_div_df <- data.frame(sample_data(pj2), alpha_diversity)
ggplot(alpha_div_df, aes(x = group, y = Shannon, fill = group)) + 
  geom_boxplot() + 
  theme_minimal() + 
  labs(title = "Shannon Diversity by Group") +
  facet_wrap(~location, scales = 'free')


#### faith pd ####

## pretreatment 
phylo_pre <- subset_samples(pj2, group == "Pre-ICI")
phylo_dist_pre <- pd(t(otu_table(phylo_pre)), phy_tree(phylo_pre),
                 include.root=F)
sample_data(phylo_pre)$PD <- phylo_dist_pre$PD

plot.pd_pre <- ggplot(sample_data(phylo_pre), aes(location, PD)) + 
  geom_boxplot() +
  xlab("Subject ID") +
  ylab("Phylogenetic Diversity")

## compare within one organ before and after treatment
phylo_dist <- pd(t(otu_table(pj2)), phy_tree(pj2),
                     include.root=F)
sample_data(pj2)$PD <- phylo_dist$PD

plot.pd_location <- ggplot(sample_data(pj2), aes(group, PD)) + 
  geom_boxplot() +
  xlab("Subject ID") +
  ylab("Phylogenetic Diversity") +
  facet_wrap(~ location)


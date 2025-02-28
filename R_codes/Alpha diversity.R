library(phyloseq)
library(picante)
library(tidyverse)

#### shannon diversity ####
load("pj2.RData")
ls()

alpha_diversity <- estimate_richness(pj2, measures = c("Shannon", "Simpson", "Chao1", "Observed"))
alpha_div_df <- data.frame(sample_data(pj2), alpha_diversity)
ggplot(alpha_div_df, aes(x = group, y = Shannon, fill = group)) + 
  geom_boxplot() + 
  theme_minimal() + 
  labs(title = "Shannon Diversity by Group") +
  facet_wrap(~location, scales = 'free')

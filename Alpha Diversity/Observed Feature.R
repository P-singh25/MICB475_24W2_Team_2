library(phyloseq)
library(picante)
library(tidyverse)
install.packages("ggsignif")
library(ggsignif)

install.packages("ggsci")  # If you haven't installed it
library(ggsci)
library(scales)

#### load phyloseq object ####
load("pj2.RData")
ls()
sample_data(pj2)
pj2_no <- subset_samples(pj2, day !="Day 0")
pj2_mod <- subset_samples(pj2_no, location !="Stool")
sample_data(pj2_mod)

#### Observed Features ####
observed_diversity <- estimate_richness(pj2_mod, measures = c("Observed"))
observed_df <- data.frame(sample_data(pj2_mod), observed_diversity)
observed_df
ggplot(observed_df, aes(x = group, y = Observed, fill = group)) + 
  geom_boxplot() + 
  theme_minimal() + 
  labs(title = "Observed Features by Group") +
  facet_wrap(~location, scales = 'free')

##Observed: Kruskall-Wallis##
kruskal.test(Observed ~ group, data = observed_df)

##Observed: ANOVA##
observed_phylo_location <- unique(observed_df$location)

for (loc in observed_phylo_location) {
  cat("\nANOVA for Location:", loc, "\n")
  subset_data <- subset(observed_df, location == loc)
  lm_pd_vs_group_log <- lm(log(Observed) ~ group, data=subset_data)
  anova_pd_vs_group_log <- aov(lm_pd_vs_group_log)
  summary(anova_pd_vs_group_log)
  print(TukeyHSD(anova_pd_vs_group_log))
}

#Observed: Significant groups
#Spleen
#Pre-ICI vs Post-ICI1
#Pre-ICI vs Post-ICI2
#Pre-ICI vs Post-ICI3
#Tumor
#Pre-ICI vs Post-ICI3
#Post-ICI1 vs Post-ICI3


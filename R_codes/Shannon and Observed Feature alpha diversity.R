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

#### shannon diversity ####
shannon_diversity <- estimate_richness(pj2_mod, measures = c("Shannon")) #, "Simpson", "Chao1", "Observed"
shannon_df <- data.frame(sample_data(pj2_mod), shannon_diversity)
ggplot(shannon_df, aes(x = group, y = Shannon, fill = group)) + 
  geom_boxplot() + 
  theme_minimal() + 
  labs(title = "Shannon Diversity by Group") +
  facet_wrap(~location, scales = 'free')

##Shannon: Kruskall-Wallis##
kruskal.test(Shannon ~ group, data = shannon_df)
shannon_df

##Shannon: ANOVA##
shannon_df_0 <- subset(shannon_df, Shannon > 0)
shannon_phylo_location <- unique(shannon_df_0$location)

for (loc in shannon_phylo_location) {
  cat("\nANOVA for Location:", loc, "\n")
  subset_data <- subset(shannon_df_0, location == loc)
  lm_pd_vs_group_log <- lm(log(Shannon) ~ group, data=subset_data)
  anova_pd_vs_group_log <- aov(lm_pd_vs_group_log)
  summary(anova_pd_vs_group_log)
  print(TukeyHSD(anova_pd_vs_group_log))
}
#Shannon: Significant groups
  #Spleen
    #Pre-ICI vs Post-ICI2
    #Pre-ICI vs Post-ICI3
  #Tumor
    #Post-ICI1 vs Post-ICI2

##Shanoon: Pre-Treatment ###
pj2_mod_pre <- subset_samples(pj2_mod, group == "Pre-ICI")
sample_data(pj2_mod_pre)

shannon_diversity_pre <- estimate_richness(pj2_mod_pre, measures = c("Shannon")) #, "Simpson", "Chao1", "Observed"
shannon_df_pre <- data.frame(sample_data(pj2_mod_pre), shannon_diversity_pre)

ggplot(shannon_df_pre, aes(x = location, y = Shannon, fill = location)) + 
  geom_boxplot(alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +  # Adds individual points
  theme_minimal() + 
  labs(title = "Shannon Diversity Across Locations", x = "Location", y = "Shannon Index") +
  theme(
    legend.position = "none",  
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels
  ) 

##Shannon: Pre-Treatment Kruskall-Wallis ##
kruskal.test(Shannon ~ location, data = shannon_df_pre)

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


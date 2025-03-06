library(phyloseq)
library(picante)
library(tidyverse)
library(ggsignif)

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
phylo_pre <- subset_samples(pj2, group == "Pre-ICI") %>%
  subset_samples(location != "Stool")

phylo_dist_pre <- pd(t(otu_table(phylo_pre)), phy_tree(phylo_pre),
                 include.root=F)
sample_data(phylo_pre)$PD <- phylo_dist_pre$PD

## Significance 
phylo_pre_sd <- sample_data(phylo_pre)
phylo_pre_wdiv <- data.frame(phylo_pre_sd, phylo_dist_pre)

### Kruskall-Wallis 
kruskal_pre <- kruskal.test(PD ~ location, data = phylo_pre_wdiv)

### ANOVA 
lm_pd_vs_loca_log <- lm(log(PD) ~ location, data=phylo_pre_wdiv)
anova_pd_vs_loca_log <- aov(lm_pd_vs_loca_log)
summary(anova_pd_vs_loca_log)
TukeyHSD(anova_pd_vs_loca_log)
### significant: Spleen-MLN, TDLN-Spleen, Tumor-Spleen

plot.pd_pre <- ggplot(sample_data(phylo_pre), aes(location, PD)) + 
  geom_boxplot() +
  geom_point() +
  xlab("Location") +
  ylab("Phylogenetic Diversity") +
  theme_classic()+
  ylim(0, 8.5) +
  geom_signif(comparisons = list(c("Spleen","MLN"), c("TDLN", "Spleen"), c("Tumor","Spleen")),
              y_position = c(8, 7, 6),
              annotations = c("0.0006","0.002","0.0004"))


## compare within one organ before and after treatment
phylo_both <- subset_samples(pj2, location != "Stool")


phylo_dist <- pd(t(otu_table(phylo_both)), phy_tree(phylo_both),
                     include.root=F)
sample_data(phylo_both)$PD <- phylo_dist$PD

## Significance 
phylo_sd <- sample_data(phylo_both)
phylo_wdiv <- data.frame(phylo_sd, phylo_dist)

### Kruskall-Wallis 
kruskal <- kruskal.test(PD ~ group, data = phylo_wdiv)

### ANOVA 

phylo_location <- unique(phylo_wdiv$location)

for (loc in phylo_location) {
  cat("\nANOVA for Location:", loc, "\n")
  
  subset_data <- subset(phylo_wdiv, location == loc)
  
  lm_pd_vs_group_log <- lm(log(PD) ~ group, data=subset_data)
  anova_pd_vs_group_log <- aov(lm_pd_vs_group_log)
  summary(anova_pd_vs_group_log)
  print(TukeyHSD(anova_pd_vs_group_log))
}

### significant groups: 
## Spleen: Pre-ICI-Post-ICI2
## MLN: Post ICI 1/2/3 - Day 0
## tumor: nothing 
## TLDN: nothing 

plot.pd_location <- ggplot(sample_data(phylo_both), aes(group, PD)) + 
  geom_boxplot() +
  xlab("Subject ID") +
  ylab("Phylogenetic Diversity") +
  facet_wrap(~ location)


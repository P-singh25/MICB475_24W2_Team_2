library(phyloseq)
library(picante)
library(tidyverse)
library(ggsignif)

install.packages("ggsci")  # If you haven't installed it
library(ggsci)

library(scales)

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

sample_data(phylo_pre)$location <- factor(sample_data(phylo_pre)$location, 
                                          levels = c("Spleen", "MLN", "TDLN", "Tumor"))  

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

plot.pd_pre <- ggplot(sample_data(phylo_pre), aes(x = location, y = PD)) + 
  geom_boxplot(aes(fill = location)) +
  geom_point() +
  xlab("Location") +
  ylab("Phylogenetic Diversity") +
  theme_classic()+
  geom_signif(comparisons = list(c("Tumor","Spleen"), c("TDLN", "Spleen"), c("Spleen","MLN")),
              y_position = c(9.5, 8.5, 7.5),
              annotations = c("P = 0.0006","P = 0.002","P = 0.0004")) +
  scale_fill_npg() +
  expand_limits(y = 10) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), 
                     labels = c("0", "2", "4", "6", "8", "10")) +
  theme(axis.text.y = element_text(size = 10, angle = 0, color = "black"), 
        axis.text.x = element_text(size = 10, angle = 0, color = "black"))

ggsave(("Alpha_Diversity_Pretreatment.png"), plot.pd_pre, width = 6, height = 3.5)



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


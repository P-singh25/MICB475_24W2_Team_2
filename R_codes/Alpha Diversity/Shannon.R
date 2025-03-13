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


##############
#Single Organ across Treatment Timeline
##############
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
#Pre-ICI vs Post-ICI2, p<0.007
#Pre-ICI vs Post-ICI3, p<0.008
#Tumor
#Post-ICI1 vs Post-ICI3, p<0.01

####Spleen
shannon_df_0_spleen <- subset(shannon_df_0, location == "Spleen")
shannon_df_0_spleen$group <- factor(shannon_df_0_spleen$group, 
                                    levels = c("Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3"))
plot.pd_spleen_Shannon <- ggplot(sample_data(shannon_df_0_spleen), aes(x = group, y = Shannon)) + 
  geom_boxplot(aes(fill = group)) +
  geom_point() +
  xlab("Group") +
  ylab("Shannon Index") +
  theme_classic()+
  geom_signif(comparisons = list(c("Pre-ICI","Post-ICI2"), c("Pre-ICI", "Post-ICI3")),
              y_position = c(4.1, 3.7),
              annotations = c("p<0.007","p<0.008")) +
  scale_fill_npg(name = "group") +
  expand_limits(y = 4.5) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), 
                     labels = c("0", "2", "4", "6", "8", "10")) +
  theme(axis.text.y = element_text(size = 10, angle = 0, color = "black"), 
        axis.text.x = element_text(size = 10, angle = 0, color = "black")) 
plot.pd_spleen_Shannon

ggsave(("Alpha_Diversity_Shannon_Spleen.png"), plot.pd_spleen_Shannon, width = 6, height = 3.5)

####Tumor
shannon_df_0_tumor <- subset(shannon_df_0, location == "Tumor")
shannon_df_0_tumor$group <- factor(shannon_df_0_tumor$group, 
                                    levels = c("Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3"))
plot.pd_tumor_Shannon <- ggplot(sample_data(shannon_df_0_tumor), aes(x = group, y = Shannon)) + 
  geom_boxplot(aes(fill = group)) +
  geom_point() +
  xlab("Group") +
  ylab("Shannon Index") +
  theme_classic()+
  geom_signif(comparisons = list(c("Post-ICI1","Post-ICI3")),
              y_position = c(4),
              annotations = c("p<0.01")) +
  scale_fill_npg(name = "group") +
  expand_limits(y = 4.5) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), 
                     labels = c("0", "2", "4", "6", "8", "10")) +
  theme(axis.text.y = element_text(size = 10, angle = 0, color = "black"), 
        axis.text.x = element_text(size = 10, angle = 0, color = "black")) 
plot.pd_tumor_Shannon
ggsave(("Alpha_Diversity_Shannon_Tumor.png"), plot.pd_tumor_Shannon, width = 6, height = 3.5)

####TDLN

shannon_df_0_TDLN <- subset(shannon_df_0, location == "TDLN")
shannon_df_0_TDLN$group <- factor(shannon_df_0_TDLN$group, 
                                   levels = c("Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3"))
plot.pd_TDLN_Shannon <- ggplot(sample_data(shannon_df_0_TDLN), aes(x = group, y = Shannon)) + 
  geom_boxplot(aes(fill = group)) +
  geom_point() +
  xlab("Group") +
  ylab("Shannon Index") +
  theme_classic()+
  scale_fill_npg(name = "group") +
  expand_limits(y = 4.5) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), 
                     labels = c("0", "2", "4", "6", "8", "10")) +
  theme(axis.text.y = element_text(size = 10, angle = 0, color = "black"), 
        axis.text.x = element_text(size = 10, angle = 0, color = "black")) 
plot.pd_TDLN_Shannon
ggsave(("Alpha_Diversity_Shannon_TDLN.png"), plot.pd_TDLN_Shannon, width = 6, height = 3.5)

####MLN
shannon_df_0_MLN <- subset(shannon_df_0, location == "MLN")
shannon_df_0_MLN$group <- factor(shannon_df_0_MLN$group, 
                                  levels = c("Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3"))
plot.pd_MLN_Shannon <- ggplot(sample_data(shannon_df_0_MLN), aes(x = group, y = Shannon)) + 
  geom_boxplot(aes(fill = group)) +
  geom_point() +
  xlab("Group") +
  ylab("Shannon Index") +
  theme_classic()+
  scale_fill_npg(name = "group") +
  expand_limits(y = 4.5) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), 
                     labels = c("0", "2", "4", "6", "8", "10")) +
  theme(axis.text.y = element_text(size = 10, angle = 0, color = "black"), 
        axis.text.x = element_text(size = 10, angle = 0, color = "black")) 
plot.pd_MLN_Shannon
ggsave(("Alpha_Diversity_Shannon_MLN.png"), plot.pd_MLN_Shannon, width = 6, height = 3.5)

###########
# Single Treatment Timepoint across Organ #
###########
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

##Shannon: Pre-Treatment ANOVA ##
shannon_df_0_pre <- subset(shannon_df_pre, Shannon > 0)
anova_test <- aov(log(Shannon) ~ location, data = shannon_df_0_pre)
summary(anova_test)
tukey_test <- TukeyHSD(anova_test)
print(tukey_test)

#Spleen_MLN, p<0.000009
#TDLN_Spleen, p<0.00005
#Tumor_Spleen, p<0.0002

#Plot with p value##
shannon_df_0_pre$location <- factor(shannon_df_0_pre$location, levels = c("Spleen", "Tumor", "TDLN", "MLN"))
plot.pd_pre_Shannon <- ggplot(sample_data(shannon_df_pre), aes(x = location, y = Shannon)) + 
  geom_boxplot(aes(fill = location)) +
  geom_point() +
  xlab("Location") +
  ylab("Shannon Index") +
  theme_classic()+
  geom_signif(comparisons = list(c("Tumor","Spleen"), c("TDLN", "Spleen"), c("Spleen","MLN")),
              y_position = c(4.5, 4.1, 3.7),
              annotations = c("p<0.0002","p<0.00005","p<0.000009")) +
  scale_fill_npg(name = "Location") +
  expand_limits(y = 5) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), 
                     labels = c("0", "2", "4", "6", "8", "10")) +
  theme(axis.text.y = element_text(size = 10, angle = 0, color = "black"), 
        axis.text.x = element_text(size = 10, angle = 0, color = "black")) 
plot.pd_pre_Shannon

ggsave(("Alpha_Diversity_Pretreatment_Shannon.png"), plot.pd_pre_Shannon, width = 6, height = 3.5)

############
##Shannon: Post-ICI1 ##
pj2_mod_post1 <- subset_samples(pj2_mod, group == "Post-ICI1")
sample_data(pj2_mod_post1)

shannon_diversity_post1 <- estimate_richness(pj2_mod_post1, measures = c("Shannon")) #, "Simpson", "Chao1", "Observed"
shannon_df_post1 <- data.frame(sample_data(pj2_mod_post1), shannon_diversity_post1)

##Shannon: Post_ICI1 Kruskall-Wallis ##
kruskal.test(Shannon ~ location, data = shannon_df_post1)
#Not Significant

##Shannon: Post_ICI1 ANOVA ##
shannon_df_0_post1 <- subset(shannon_df_post1, Shannon > 0)
anova_test_post1 <- aov(log(Shannon) ~ location, data = shannon_df_0_post1)
summary(anova_test_post1)
tukey_test_post1 <- TukeyHSD(anova_test_post1)
print(tukey_test_post1)

#Spleen_MLN, p<0.55
#TDLN_Spleen, p<0.54
#Tumor_Spleen, p<0.22

#Plot with p value##
shannon_df_0_post1$location <- factor(shannon_df_0_post1$location, levels = c("Spleen", "Tumor", "TDLN", "MLN"))
plot.pd_post1_Shannon <- ggplot(sample_data(shannon_df_post1), aes(x = location, y = Shannon)) + 
  geom_boxplot(aes(fill = location)) +
  geom_point() +
  xlab("Location") +
  ylab("Shannon Index") +
  theme_classic()+
  geom_signif(comparisons = list(c("Tumor","Spleen"), c("TDLN", "Spleen"), c("Spleen","MLN")),
              y_position = c(4.8, 4.3, 3.9),
              annotations = c("p<0.22","p<0.54","p<0.55")) +
  scale_fill_npg(name = "Location") +
  expand_limits(y = 5.3) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), 
                     labels = c("0", "2", "4", "6", "8", "10")) +
  theme(axis.text.y = element_text(size = 10, angle = 0, color = "black"), 
        axis.text.x = element_text(size = 10, angle = 0, color = "black")) 
plot.pd_post1_Shannon

ggsave(("Alpha_Diversity_Pretreatment_Shannon_Post-ICI1.png"), plot.pd_post1_Shannon, width = 6, height = 3.5)

#############
##Shannon: Post-ICI2 ##
pj2_mod_post2 <- subset_samples(pj2_mod, group == "Post-ICI2")
sample_data(pj2_mod_post2)

shannon_diversity_post2 <- estimate_richness(pj2_mod_post2, measures = c("Shannon")) #, "Simpson", "Chao1", "Observed"
shannon_df_post2 <- data.frame(sample_data(pj2_mod_post2), shannon_diversity_post2)

##Shannon: Post_ICI2 Kruskall-Wallis ##
kruskal.test(Shannon ~ location, data = shannon_df_post2)
#Not Significant

##Shannon: Post_ICI2 ANOVA ##
shannon_df_0_post2 <- subset(shannon_df_post2, Shannon > 0)
anova_test_post2 <- aov(log(Shannon) ~ location, data = shannon_df_0_post2)
summary(anova_test_post2)
tukey_test_post2 <- TukeyHSD(anova_test_post2)
print(tukey_test_post2)

#Spleen_MLN, p<0.26
#TDLN_Spleen, p<0.44
#Tumor_Spleen, p<0.18

#Plot with p value##
shannon_df_0_post2$location <- factor(shannon_df_0_post2$location, levels = c("Spleen", "Tumor", "TDLN", "MLN"))
plot.pd_post2_Shannon <- ggplot(sample_data(shannon_df_0_post2), aes(x = location, y = Shannon)) + 
  geom_boxplot(aes(fill = location)) +
  geom_point() +
  xlab("Location") +
  ylab("Shannon Index") +
  theme_classic()+
  geom_signif(comparisons = list(c("Tumor","Spleen"), c("TDLN", "Spleen"), c("Spleen","MLN")),
              y_position = c(4.8, 4.3, 3.9),
              annotations = c("p<0.18","p<0.44","p<0.26")) +
  scale_fill_npg(name = "Location") +
  expand_limits(y = 5.3) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), 
                     labels = c("0", "2", "4", "6", "8", "10")) +
  theme(axis.text.y = element_text(size = 10, angle = 0, color = "black"), 
        axis.text.x = element_text(size = 10, angle = 0, color = "black")) 
plot.pd_post2_Shannon

ggsave(("Alpha_Diversity_Pretreatment_Shannon_Post-ICI2.png"), plot.pd_post2_Shannon, width = 6, height = 3.5)



###############
##Shannon: Post-ICI3 ##
pj2_mod_post <- subset_samples(pj2_mod, group == "Post-ICI3")
sample_data(pj2_mod_post)

shannon_diversity_post <- estimate_richness(pj2_mod_post, measures = c("Shannon")) #, "Simpson", "Chao1", "Observed"
shannon_df_post <- data.frame(sample_data(pj2_mod_post), shannon_diversity_post)

##Shannon: Post_ICI3 Kruskall-Wallis ##
kruskal.test(Shannon ~ location, data = shannon_df_post)
#Not Significant

##Shannon: Post_ICI3 ANOVA ##
shannon_df_0_post <- subset(shannon_df_0_post, Shannon > 0)
anova_test_post <- aov(log(Shannon) ~ location, data = shannon_df_0_post)
summary(anova_test_post)
tukey_test_post <- TukeyHSD(anova_test_post)
print(tukey_test_post)

#Spleen_MLN, p<0.96
#TDLN_Spleen, p<0.62
#Tumor_Spleen, p<0.20

#Plot with p value##
shannon_df_0_post$location <- factor(shannon_df_0_post$location, levels = c("Spleen", "Tumor", "TDLN", "MLN"))
plot.pd_post_Shannon <- ggplot(sample_data(shannon_df_0_post), aes(x = location, y = Shannon)) + 
  geom_boxplot(aes(fill = location)) +
  geom_point() +
  xlab("Location") +
  ylab("Shannon Index") +
  theme_classic()+
  geom_signif(comparisons = list(c("Tumor","Spleen"), c("TDLN", "Spleen"), c("Spleen","MLN")),
              y_position = c(4.8, 4.3, 3.9),
              annotations = c("p<0.20","p<0.62","p<0.96")) +
  scale_fill_npg(name = "Location") +
  expand_limits(y = 5.3) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), 
                     labels = c("0", "2", "4", "6", "8", "10")) +
  theme(axis.text.y = element_text(size = 10, angle = 0, color = "black"), 
        axis.text.x = element_text(size = 10, angle = 0, color = "black")) 
plot.pd_post_Shannon

ggsave(("Alpha_Diversity_Pretreatment_Shannon_Post-ICI3.png"), plot.pd_post_Shannon, width = 6, height = 3.5)




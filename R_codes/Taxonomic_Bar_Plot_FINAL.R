# Load required packages
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyverse)



load("pj2.RData")
pj2_melted <- psmelt(pj2)

#### test on spleen only ####
pj2_melted_spleen <- pj2_melted %>%
  filter(location == "Spleen") %>%
  filter(group != "Day0") %>%
  mutate(group = factor(group, levels = c("Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3"))) %>%
  drop_na()

pj2_group_total <- pj2_melted_spleen %>%
  group_by(group) %>%
  mutate(Total_Abundance = sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

pj2_ra <- pj2_group_total %>%
  group_by(group, OTU) %>%
  mutate(Relative_Abundance = sum(Abundance)/Total_Abundance)

taxa_spleen <- ggplot(data = pj2_ra, aes(x = group, y = Relative_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
        axis.text.y = element_text(size = 10, angle = 0, color = "black"), 
        axis.text.x = element_text(size = 10, angle = 0, color = "black"), 
        axis.line = element_line(color = "black", size = 0), 
        axis.title.x = element_text(size = 12, face = "bold"),  # Bigger x-axis label
        axis.title.y = element_text(size = 12, face = "bold"), 
        legend.title = element_text(size = 10, face = "bold", hjust = 0.5))+
  guides(fill = guide_legend(ncol = 1)) +
  scale_fill_viridis_d(option = "turbo") +
  labs(x = "Treatment", y = "Relative Abundance", fill = "Bacterial Phylum")
  

ggsave("Spleen_Taxa_Bar.png", taxa_spleen, width = 7, height = 10)


#### full #### 

## melt the phyloseq object to data frame
pj2_melted <- pj2_melted %>%
  filter(group != "Day0") %>% 
  mutate(group = factor(group, levels = c("Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3"))) %>%
  drop_na()

## calculate the total number of reads in each group of each location 
pj2_all_total <- pj2_melted %>%
  group_by(location, group) %>%
  mutate(Total_Abundance = sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

## calculate relative abundance 
pj2_all_ra <- pj2_all_total %>%
  group_by(group, OTU, location) %>%
  mutate(Relative_Abundance = sum(Abundance)/Total_Abundance)

taxa_plot <- ggplot(data = pj2_all_ra, aes(x = group, y = Relative_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~factor(location), scales = "free_x", nrow = 1) +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(size = 12, face = "bold"),  # Bigger x-axis label
        axis.title.y = element_text(size = 12, face = "bold"), 
        legend.title = element_text(size = 10, face = "bold", hjust = 0.5), 
        strip.text = element_text(size = 14, face = "bold"))+
  guides(fill = guide_legend(ncol = 1)) +
  scale_fill_viridis_d(option = "turbo") +
  labs(x = "Treatment", y = "Relative Abundance (%)", fill = "Bacterial Phylum")

ggsave("Taxa_Bar.png", taxa_plot, width = 15, height = 10)




#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load data ####
load("pj2.RData")

pj2_spleen <- pj2 %>%
  subset_samples(location == "Spleen") %>%
  subset_samples(group != "Day 0")

#### Indicator Species/Taxa Analysis ####
# glom to Genus
pj2_glom <- tax_glom(pj2_spleen, "Genus", NArm = FALSE)
pj2_glom_RA <- transform_sample_counts(pj2_glom, fun=function(x) x/sum(x))

sum(is.na(otu_table(pj2_glom_RA)))  # Count NA values
otu_table(pj2_glom_RA)[is.na(otu_table(pj2_glom_RA))] <- 0

sum(is.na(sample_data(pj2_glom_RA)$group))  # Count NA values in group

#ISA
pj2_is <- multipatt(t(otu_table(pj2_glom_RA)), cluster = sample_data(pj2_glom_RA)$`group`)
summary(pj2_is)
taxtable <- tax_table(pj2) %>% as.data.frame() %>% rownames_to_column(var="ASV")

## filter for significant results 
isa_results_spleen <- pj2_is$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05)

## save the table as csv 
write.csv(isa_results_spleen, "indicator_species_spleen.csv", row.names = FALSE)

#### filter further based on stat 
isa_results_spleen_stri <- isa_results_spleen %>%
  filter(stat > 0.7) %>%
  as.data.frame()

## save the table as csv 
write.csv(isa_results_spleen_stri, "indicator_species_spleen_stringent.csv", row.names = FALSE)

## remove the columns not needed for data presentation
df_filtered <- isa_results_spleen_stri[, !names(isa_results_spleen_stri) %in% c("ASV", "index", "Species")]

## save the table as csv                                    
write.csv(df_filtered, "indicator_species_spleen_filter.csv", row.names = FALSE)


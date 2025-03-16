#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load data ####
load("pj2.RData")

#### Indicator Species/Taxa Analysis ####
# glom to Genus
pj2_glom <- tax_glom(pj2, "Genus", NArm = FALSE)
pj2_glom_RA <- transform_sample_counts(pj2_glom, fun=function(x) x/sum(x))

sum(is.na(otu_table(pj2_glom_RA)))  # Count NA values
otu_table(pj2_glom_RA)[is.na(otu_table(pj2_glom_RA))] <- 0

sum(is.na(sample_data(pj2_glom_RA)$group))  # Count NA values in group

#ISA
pj2_is <- multipatt(t(otu_table(pj2_glom_RA)), cluster = sample_data(pj2_glom_RA)$`group`)
summary(pj2_is)
taxtable <- tax_table(pj2) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# consider that your table is only going to be resolved up to the genus level, be wary of 
# anything beyond the glomed taxa level
isa_results <- pj2_is$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% View()

write.csv(isa_results, "indicator_species_results.csv", row.names = TRUE)

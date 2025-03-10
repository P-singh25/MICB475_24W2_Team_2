#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(DESeq2)

#### Load data ####
load("pj2.RData")

#### DESeq for spleen ####
pj2_spleen <- pj2 %>%
  subset_samples(location == "Spleen")

pj2_deseq_spleen <- phyloseq_to_deseq2(pj2_spleen, ~`group`)
DESEQ_pj2_spleen <- DESeq(pj2_deseq_spleen)

res <- results(DESEQ_pj2_spleen, contrast = c("group", "Post-ICI1", "Pre-ICI"))

View(res)

# Look at results 

## Volcano plot: effect size VS significance
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

# Remove rows with NA in padj
res_clean <- res[!is.na(res$padj), ]

# Now, filter for upregulated species
upregulated_species <- res_clean[res_clean$log2FoldChange > 0 & res_clean$padj < 0.05, ]

# View the upregulated species
View(upregulated_species)

# Sort by highest fold change
upregulated_species <- upregulated_species[order(-upregulated_species$log2FoldChange), ]

# Print top 10 most upregulated taxa
head(upregulated_species, 10)
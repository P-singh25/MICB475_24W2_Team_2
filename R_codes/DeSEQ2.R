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

res <- results(DESEQ_pj2_spleen, tidy=TRUE, contrast = c("group", "Post-ICI1", "Pre-ICI"))

# To get table of results
sigASVs <- res %>% 
  filter(padj<0.05 & log2FoldChange>2) %>%
  dplyr::rename(ASV=row)

# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

# Prune phyloseq file
mpt_DESeq <- prune_taxa(sigASVs_vec,pj2_spleen)

sigASVs <- tax_table(mpt_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

deseq_spleen_plot <- ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggsave("Deseq_Spleen.png", deseq_spleen_plot, width = 6, height = 3.5)

#### all locations ####
pj2_all <- pj2 %>%
  subset_samples(location != "Stool")

## zero errors occur 
## pj2_deseq_all <- phyloseq_to_deseq2(pj2_all, ~`group`)
## DESEQ_pj2_all <- DESeq(pj2_deseq_all)

pj2_plus1 <- transform_sample_counts(pj2_all, function(x) x+1)
pj2_all_deseq <- phyloseq_to_deseq2(pj2_plus1, ~`group`)
DESEQ_mpt <- DESeq(pj2_all_deseq)
res_all <- results(DESEQ_mpt, tidy=TRUE, 
               #this will ensure that No is your reference group
               contrast = c("group", "Post-ICI1", "Pre-ICI")) 
res_all <- res_all[!is.na(res_all$padj), ]

# To get table of results
sigASVs_all <- res_all %>% 
  filter(padj<0.05 & log2FoldChange>4) %>%
  filter(!is.na(padj))%>%
  dplyr::rename(ASV=row)

# Get only asv names
sigASVs_vec_all <- sigASVs_all %>%
  pull(ASV)

# Prune phyloseq file
pj2_all_DESeq <- prune_taxa(sigASVs_vec_all,pj2_all)

sigASVs_all <- tax_table(pj2_all_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_all) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  filter(!is.na(Genus) & Genus != "NA.1" & Genus != "NA.2")

deseq_all_plot <- ggplot(sigASVs_all) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggsave("Deseq_all.png", deseq_all_plot, width = 6, height = 3.5)



#### from Ran #### 
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

## Volcano plot: effect size VS significance
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))


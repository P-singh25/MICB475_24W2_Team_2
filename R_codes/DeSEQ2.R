#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(DESeq2)

#### Load data ####
load("pj2.RData")

#### DESeq for spleen ####
pj2_spleen <- pj2 %>%
  subset_samples(location == "Spleen") %>%
  subset_samples(group != "Day0")

pj2_deseq_spleen <- phyloseq_to_deseq2(pj2_spleen, ~`group`)
DESEQ_pj2_spleen <- DESeq(pj2_deseq_spleen)
res_spleen_post1 <- results(DESEQ_pj2_spleen, tidy=TRUE, 
                            contrast = c("group", "Post-ICI1", "Pre-ICI"))
View(res_spleen_post1)

# To get table of results
sigASVs_spleen_post1 <- res_spleen_post1 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get only asv names
sigASVs_name_spleen_post1 <- sigASVs_spleen_post1 %>%
  pull(ASV)

# Prune phyloseq file
spleen_post1_DESeq <- prune_taxa(sigASVs_name_spleen_post1, pj2_spleen)

sigASVs_spleen_post1 <- tax_table(spleen_post1_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_spleen_post1) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

deseq_spleen_post1_plot <- ggplot(sigASVs_spleen_post1) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = Phylum), stat="identity")+
  # geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), size = 0.5, width = 0.5)  +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1, color = "black"), 
        axis.text.y = element_text(color = "black"))

ggsave("Deseq_Spleen_prevspost1.png", deseq_spleen_post1_plot, width = 6, height = 3.5)

#### Spleen post 2 ####
res_spleen_post2 <- results(DESEQ_pj2_spleen, tidy=TRUE, 
                            contrast = c("group", "Post-ICI2", "Pre-ICI"))

# To get table of results
sigASVs_spleen_post2 <- res_spleen_post2 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get only asv names
sigASVs_name_spleen_post2 <- sigASVs_spleen_post2 %>%
  pull(ASV)

# Prune phyloseq file
spleen_post2_DESeq <- prune_taxa(sigASVs_name_spleen_post2, pj2_spleen)

sigASVs_spleen_post2 <- tax_table(spleen_post2_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_spleen_post2) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

deseq_spleen_post2_plot <- ggplot(sigASVs_spleen_post2) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = Phylum), stat="identity")+
  # geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), size = 0.5, width = 0.5)  +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1, color = "black"), 
        axis.text.y = element_text(color = "black"))

ggsave("Deseq_Spleen_prevspost2.png", deseq_spleen_post2_plot, width = 6, height = 3.5)


#### spleen post 3 ####
res_spleen_post3 <- results(DESEQ_pj2_spleen, tidy=TRUE, 
                            contrast = c("group", "Post-ICI3", "Pre-ICI"))

# To get table of results
sigASVs_spleen_post3 <- res_spleen_post3 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get only asv names
sigASVs_name_spleen_post3 <- sigASVs_spleen_post3 %>%
  pull(ASV)

# Prune phyloseq file
spleen_post3_DESeq <- prune_taxa(sigASVs_name_spleen_post3, pj2_spleen)

sigASVs_spleen_post3 <- tax_table(spleen_post3_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_spleen_post3) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

deseq_spleen_post3_plot <- ggplot(sigASVs_spleen_post3) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = Phylum), stat="identity")+
  # geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), size = 0.5, width = 0.5)  +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1, color = "black"), 
        axis.text.y = element_text(color = "black"))

ggsave("Deseq_Spleen_prevspost3.png", deseq_spleen_post3_plot, width = 6, height = 3.5)





#### tumor post 1 ####
pj2_tumor <- pj2 %>%
  subset_samples(location == "Tumor") %>%
  subset_samples(group != "Day0")

pj2_deseq_tumor <- phyloseq_to_deseq2(pj2_tumor, ~`group`)
DESEQ_pj2_tumor <- DESeq(pj2_deseq_tumor)
res_tumor_post1 <- results(DESEQ_pj2_tumor, tidy=TRUE, 
                            contrast = c("group", "Post-ICI1", "Pre-ICI"))

# To get table of results
sigASVs_tumor_post1 <- res_tumor_post1 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get only asv names
sigASVs_name_tumor_post1 <- sigASVs_tumor_post1 %>%
  pull(ASV)

# Prune phyloseq file
tumor_post1_DESeq <- prune_taxa(sigASVs_name_tumor_post1, pj2_tumor)

sigASVs_tumor_post1 <- tax_table(tumor_post1_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_tumor_post1) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  filter(Genus != "NA")

deseq_tumor_post1_plot <- ggplot(sigASVs_tumor_post1) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = Phylum), stat="identity")+
  # geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), size = 0.5, width = 0.5)  +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1, color = "black"), 
        axis.text.y = element_text(color = "black"))

ggsave("Deseq_Tumor_prevspost1.png", deseq_tumor_post1_plot, width = 6, height = 3.5)





#### TDLN post 1 ####
pj2_TDLN <- pj2 %>%
  subset_samples(location == "TDLN") %>%
  subset_samples(group != "Day0") %>%
  transform_sample_counts(function(x) x+1)

pj2_deseq_TDLN <- phyloseq_to_deseq2(pj2_TDLN, ~`group`)
DESEQ_pj2_TDLN <- DESeq(pj2_deseq_TDLN)
res_TDLN_post1 <- results(DESEQ_pj2_TDLN, tidy=TRUE, 
                            contrast = c("group", "Post-ICI1", "Pre-ICI"))

# To get table of results
sigASVs_TDLN_post1 <- res_TDLN_post1 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get only asv names
sigASVs_name_TDLN_post1 <- sigASVs_TDLN_post1 %>%
  pull(ASV)

# Prune phyloseq file
TDLN_post1_DESeq <- prune_taxa(sigASVs_name_TDLN_post1, pj2_TDLN)

sigASVs_TDLN_post1 <- tax_table(TDLN_post1_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_TDLN_post1) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  filter(Genus != "NA")

deseq_TDLN_post1_plot <- ggplot(sigASVs_TDLN_post1) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = Phylum), stat="identity")+
  # geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), size = 0.5, width = 0.5)  +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1, color = "black"), 
        axis.text.y = element_text(color = "black"))

ggsave("Deseq_TDLN_prevspost1.png", deseq_TDLN_post1_plot, width = 6, height = 3.5)





#### MLN ####
pj2_MLN <- pj2 %>%
  subset_samples(location == "MLN") %>%
  subset_samples(group != "Day0") %>%
  transform_sample_counts(function(x) x+1)

pj2_deseq_MLN <- phyloseq_to_deseq2(pj2_MLN, ~`group`)
DESEQ_pj2_MLN <- DESeq(pj2_deseq_MLN)
res_MLN_post1 <- results(DESEQ_pj2_MLN, tidy=TRUE, 
                          contrast = c("group", "Post-ICI1", "Pre-ICI"))

# To get table of results
sigASVs_MLN_post1 <- res_MLN_post1 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get only asv names
sigASVs_name_MLN_post1 <- sigASVs_MLN_post1 %>%
  pull(ASV)

# Prune phyloseq file
MLN_post1_DESeq <- prune_taxa(sigASVs_name_MLN_post1, pj2_MLN)

sigASVs_MLN_post1 <- tax_table(MLN_post1_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_MLN_post1) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  filter(Genus != "NA")

deseq_MLN_post1_plot <- ggplot(sigASVs_MLN_post1) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = Phylum), stat="identity")+
  # geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), size = 0.5, width = 0.5)  +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1, color = "black"), 
        axis.text.y = element_text(color = "black"))

ggsave("Deseq_MLN_prevspost1.png", deseq_MLN_post1_plot, width = 6, height = 3.5)




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




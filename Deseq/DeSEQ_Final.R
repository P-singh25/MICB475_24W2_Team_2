#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(DESeq2)

#### Load data ####
load("pj2.RData")

#### DESeq for spleen pre vs post3 ####
pj2_spleen <- pj2 %>%
  subset_samples(location == "Spleen") %>%
  subset_samples(group != "Day0")

pj2_deseq_spleen <- phyloseq_to_deseq2(pj2_spleen, ~`group`)
DESEQ_pj2_spleen <- DESeq(pj2_deseq_spleen)

res_spleen_post3 <- results(DESEQ_pj2_spleen, tidy=TRUE, contrast = c("group", "Post-ICI3", "Pre-ICI"))

# To get table of results
sigASVs_spleen_post3 <- res_spleen_post3 %>% 
  filter(padj<0.05 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get only asv names
sigASVs_vec_spleen_post3 <- sigASVs_spleen_post3 %>%
  pull(ASV)

# Prune phyloseq file
spleen_post3_DESeq <- prune_taxa(sigASVs_vec_spleen_post3,pj2_spleen)

sigASVs_spleen_post3 <- tax_table(spleen_post3_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_spleen_post3) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

# visualze the data using box plot 
deseq_spleen_post3_plot <- ggplot(sigASVs_spleen_post3) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = Phylum), stat="identity", width = 0.7)+
  # geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), size = 0.5, width = 0.5)  +
  theme_bw() +
  scale_x_discrete(labels = c("g__Staphylococcus" = "Staphylococcus", 
                              "g__Phascolarctobacterium" = "Phascolarctobacterium", 
                              "g__Enterococcus" = "Enterococcus", 
                              "g__Akkermansia" = "Akkermansia", 
                              "g__Dubosiella" = "Dubosiella", 
                              "g__Clostridioides" = "Clostridioides"), 
                   expand = c(0.1, 0)) +
  scale_fill_manual(values = c("p__Firmicutes" = "#7fc680", "p__Verrucomicrobiota" = "#0e9ee0"),  # Set colors
                    labels = c("p__Firmicutes" = "Firmicutes", "p__Verrucomicrobiota" = "Verrucomicrobiota"))+
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.2) +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1, color = "black"), 
        axis.text.y = element_text(color = "black"), 
        legend.title = element_text(size = 12, hjust = 0.5))

ggsave("Deseq_Spleen_prevspost3.png", deseq_spleen_post3_plot, width = 6, height = 3.5)



#### Tumor post ICI 3 ####
pj2_tumor <- pj2 %>%
  subset_samples(location == "Tumor") %>%
  subset_samples(group != "Day0")

pj2_deseq_tumor <- phyloseq_to_deseq2(pj2_tumor, ~`group`)
DESEQ_pj2_tumor <- DESeq(pj2_deseq_tumor)

res_tumor <- results(DESEQ_pj2_tumor, tidy=TRUE, contrast = c("group", "Post-ICI3", "Pre-ICI"))

# To get table of results
sigASVs_tumor <- res_tumor %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get only asv names
sigASVs_tumor_vec <- sigASVs_tumor %>%
  pull(ASV)

# Prune phyloseq file
tumor_DESeq <- prune_taxa(sigASVs_tumor_vec,pj2_tumor)

sigASVs_tumor <- tax_table(tumor_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_tumor) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  filter(Genus != "NA")

deseq_tumor_plot <- ggplot(sigASVs_tumor) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = Phylum), stat="identity")+
  # geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), size = 0.5, width = 0.5)  +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, color = "black"), 
        axis.text.y = element_text(color = "black"))

ggsave("Deseq_Tumor_prevspost3.png", deseq_tumor_plot, width = 6, height = 3.5)


#### TDLN post ICI 3####
pj2_TDLN <- pj2 %>%
  subset_samples(location == "TDLN") %>%
  subset_samples(group != "Day0")

# transform the data 
pj2_TDLN_plus1 <- transform_sample_counts(pj2_TDLN, function(x) x+1)

pj2_deseq_TDLN <- phyloseq_to_deseq2(pj2_TDLN_plus1, ~`group`)
DESEQ_pj2_TDLN <- DESeq(pj2_deseq_TDLN)

res_TDLN <- results(DESEQ_pj2_TDLN, tidy=TRUE, 
                    contrast = c("group", "Post-ICI3", "Pre-ICI"))

# To get table of results
sigASVs_TDLN <- res_TDLN %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get only asv names
sigASVs_TDLN_vec <- sigASVs_TDLN %>%
  pull(ASV)

# Prune phyloseq file
TDLN_DESeq <- prune_taxa(sigASVs_TDLN_vec,pj2_TDLN)

sigASVs_TDLN <- tax_table(TDLN_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_TDLN) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  filter(Genus != "NA")

deseq_TDLN_plot <- ggplot(sigASVs_TDLN) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = Phylum), stat="identity")+
  # geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), size = 0.5, width = 0.5)  +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, color = "black"), 
        axis.text.y = element_text(color = "black"))

ggsave("Deseq_TDLN_prevspost3.png", deseq_TDLN_plot, width = 6, height = 3.5)


#### MLN post ICI 3 ####
pj2_MLN <- pj2 %>%
  subset_samples(location == "MLN") %>%
  subset_samples(group != "Day0")

pj2_MLN_plus1 <- transform_sample_counts(pj2_MLN, function(x) x+1)

pj2_deseq_MLN <- phyloseq_to_deseq2(pj2_MLN_plus1, ~`group`)
DESEQ_pj2_MLN <- DESeq(pj2_deseq_MLN)

res_MLN <- results(DESEQ_pj2_MLN, tidy=TRUE, 
                   contrast = c("group", "Post-ICI3", "Pre-ICI"))

# To get table of results
sigASVs_MLN <- res_MLN %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get only asv names
sigASVs_MLN_vec <- sigASVs_MLN %>%
  pull(ASV)

# Prune phyloseq file
MLN_DESeq <- prune_taxa(sigASVs_MLN_vec,pj2_MLN)

sigASVs_MLN <- tax_table(MLN_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_MLN) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  filter(Genus != "NA")

deseq_MLN_plot <- ggplot(sigASVs_MLN) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = Phylum), stat="identity")+
  # geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), size = 0.5, width = 0.5)  +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, color = "black"), 
        axis.text.y = element_text(color = "black"))

ggsave("Deseq_MLN_prevspost3.png", deseq_MLN_plot, width = 6, height = 3.5)



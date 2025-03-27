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

res <- results(DESEQ_pj2_spleen, tidy=TRUE, contrast = c("group", "Post-ICI1", "Pre-ICI"))

# To get table of results
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

# Prune phyloseq file
mpt_DESeq <- prune_taxa(sigASVs_vec,pj2_spleen)

sigASVs_spleen <- tax_table(mpt_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_spleen) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

deseq_spleen_plot <- ggplot(sigASVs_spleen) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = Phylum), stat="identity")+
  # geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), size = 0.5, width = 0.5)  +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1, color = "black"), 
        axis.text.y = element_text(color = "black"))

ggsave("Deseq_Spleen_prevspost1.png", deseq_spleen_plot, width = 6, height = 3.5)


#### DESeq for spleen pre vs post2 ####
pj2_spleen <- pj2 %>%
  subset_samples(location == "Spleen") %>%
  subset_samples(group != "Day0")

pj2_deseq_spleen <- phyloseq_to_deseq2(pj2_spleen, ~`group`)
DESEQ_pj2_spleen <- DESeq(pj2_deseq_spleen)

res_spleen_post2 <- results(DESEQ_pj2_spleen, tidy=TRUE, contrast = c("group", "Post-ICI2", "Pre-ICI"))

# To get table of results
sigASVs_spleen_post2 <- res_spleen_post2 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get only asv names
sigASVs_vec_spleen_post2 <- sigASVs_spleen_post2 %>%
  pull(ASV)

# Prune phyloseq file
spleen_post2_DESeq <- prune_taxa(sigASVs_vec_spleen_post2,pj2_spleen)

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


#### DESeq for spleen pre vs post3 ####
pj2_spleen <- pj2 %>%
  subset_samples(location == "Spleen") %>%
  subset_samples(group != "Day0")

pj2_deseq_spleen <- phyloseq_to_deseq2(pj2_spleen, ~`group`)
DESEQ_pj2_spleen <- DESeq(pj2_deseq_spleen)

res_spleen_post3 <- results(DESEQ_pj2_spleen, tidy=TRUE, contrast = c("group", "Post-ICI3", "Pre-ICI"))

# To get table of results
sigASVs_spleen_post3 <- res_spleen_post3 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
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

deseq_spleen_post3_plot <- ggplot(sigASVs_spleen_post3) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = Phylum), stat="identity")+
  # geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), size = 0.5, width = 0.5)  +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1, color = "black"), 
        axis.text.y = element_text(color = "black"))

ggsave("Deseq_Spleen_prevspost3.png", deseq_spleen_post3_plot, width = 6, height = 3.5)



#### Tumor post 3 ####
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


#### TDLN post 3####
pj2_TDLN <- pj2 %>%
  subset_samples(location == "TDLN") %>%
  subset_samples(group != "Day0")

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


#### MLN post 3 ####
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

# To get table of results
sigASVs_all <- res_all %>% 
  filter(padj<0.01 & log2FoldChange>4) %>%
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


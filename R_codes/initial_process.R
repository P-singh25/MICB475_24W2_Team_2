library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

# Load data
metadata <- "pj2_export/melanoma_metadata.tsv"
meta_pre <- read_delim(metadata, delim="\t")

table <- "pj2_export/feature-table.txt"
otu <- read_delim(file = table, delim="\t", skip=1)

taxonomy <- "pj2_export/taxonomy.tsv"
tax <- read_delim(taxonomy, delim="\t")

phylotree_data <- "pj2_export/tree.nwk"
phylotree <- read.tree(phylotree_data)

#optimizing the metadata file
map_day_label <- function(day_str) {
  day_num <- as.numeric(str_remove(day_str, "Day "))  # Extract numeric part
  
  if (day_num == 0) {
    return("Day 0")
  } else if (day_num >= 1 & day_num <= 4) {
    return("Pre-ICI")
  } else if (day_num >= 5 & day_num <= 8) {
    return("Post-ICI1")
  } else if (day_num >= 9 & day_num <= 12) {
    return("Post-ICI2")
  } else if (day_num >= 13 & day_num <= 15) {
    return("Post-ICI3")
  } else {
    return(day_str)  # Keep original if unexpected value
  }
}

# Apply the mapping function to the 'day' column
meta <- meta_pre %>%
  mutate(group = map_chr(day, map_day_label))

# Save the updated file
write_tsv(meta, "pj2_export/melanoma_metadata_updated.tsv")

#Format OTU table as matrix, where rownames and colnames as OTUs and sampleIDs, respectively
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

#Format sample metadata
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$`sample-id`
SAMP <- sample_data(samp_df)

#Formatting taxonomy
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)

#Make phyloseq object
pj2 <- phyloseq(OTU, SAMP, TAX, phylotree)
save(mt, file="pj2.RData")

# Remove non-bacterial sequences, if any
pj2_filt <- subset_taxa(pj2,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# Remove ASVs that have less than 5 counts total
pj2_filt_nolow <- filter_taxa(pj2_filt, function(x) sum(x)>5, prune = TRUE)

# Remove samples with less than 100 reads
pj2_filt_nolow_samps <- prune_samples(sample_sums(pj2_filt_nolow)>100, pj2_filt_nolow)

#check maximum sequencing depth
max_depth <- max(sample_sums(pj2_filt_nolow_samps))
print(max_depth)

#Look at the rarefaction curve to choose rarefaction depth
#rarecurve(t(as.data.frame(otu_table(pj2_filt_nolow_samps))), cex=0.1)

#Rarefy with a sequencing depth of 1000
rare <- rarefy_even_depth(pj2_filt_nolow_samps, rngseed = 1, sample.size = 650) 
data_rare <- as.data.frame(sample_data(rare))


#alpha diversity
alpha_diversity <- estimate_richness(rare, measures = c("Shannon", "Simpson", "Chao1", "Observed"))
alpha_div_df <- data.frame(sample_data(rare), alpha_diversity)
ggplot(alpha_div_df, aes(x = group, y = Shannon, fill = group)) + 
  geom_boxplot() + 
  theme_minimal() + 
  labs(title = "Shannon Diversity by Group") +
  facet_wrap(~location, scales = 'free')

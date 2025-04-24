# load library
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(ggplot2)

#
sample_data(pj2)

# Subset the phyloseq object for MLN Pre-ICI and MLN Post-ICI3
Spleen_Pre <- subset_samples(pj2, location == "Spleen" & group == "Pre-ICI")
Spleen_Post1 <- subset_samples(pj2, location == "Spleen" & group == "Post-ICI1")
Spleen_Post2 <- subset_samples(pj2, location == "Spleen" & group == "Post-ICI2")
Spleen_Post3 <- subset_samples(pj2, location == "Spleen" & group == "Post-ICI3")


# Transform sample counts to relative abundance (optional but recommended for core analysis)

Spleen_Pre_relabund <- transform_sample_counts(Spleen_Pre, fun = function(x) x / sum(x))
Spleen_Post1_relabund <- transform_sample_counts(Spleen_Post1, fun = function(x) x / sum(x))
Spleen_Post2_relabund <- transform_sample_counts(Spleen_Post2, fun = function(x) x / sum(x))
Spleen_Post3_relabund <- transform_sample_counts(Spleen_Post3, fun = function(x) x / sum(x))

# Identify core microbiome for each group
# Core microbiome is defined as taxa present in at least 80% of samples (prevalence = 0.3) with any non-zero abundance (detection = 0)


core_Spleen_Pre <- core_members(Spleen_Pre_relabund, detection = 0, prevalence = 0.3)
core_Spleen_Post1 <- core_members(Spleen_Post1_relabund, detection = 0, prevalence = 0.3)
core_Spleen_Post2 <- core_members(Spleen_Post2_relabund, detection = 0, prevalence = 0.3)
core_Spleen_Post3 <- core_members(Spleen_Post3_relabund, detection = 0, prevalence = 0.3)
########## Post ICI-3 ############

# Find shared taxa 
shared_taxa_post <- Reduce(intersect, list(core_Spleen_Pre, core_Spleen_Post1, core_Spleen_Post2, core_Spleen_Post3))
# Prune phyloseq object to shared taxa
shared_physeq_post <- prune_taxa(shared_taxa_post, pj2)
# View taxonomic information
tax_table(shared_physeq_post)

# Plot the venn diagram
core_all_post_ICI3 <- ggVennDiagram(
  x = list(
    "Pre-ICI" = core_Spleen_Pre,
    "Post-ICI1" = core_Spleen_Post1,
    "Post-ICI2" = core_Spleen_Post2,
    "Post-ICI3" = core_Spleen_Post3
  ), label_alpha = 0, set_size = 6, label_size = 6
  #, label = "count"
)+
  scale_fill_gradient(low = "white", high = "orange") +  # Fill colors from white to orange
  scale_color_manual(values = c("black", "black", "black", "black")) +
  theme(text = element_text(size = 13))
 # theme(legend.position = "none")

core_all_post_ICI3

# Save the venn diagram
ggsave(("Core_microbiome_Spleen.png"), core_all_post_ICI3, width = 10, height = 10)

#list unique from post-ICI3
unique_post3 <- setdiff(core_Spleen_Post3, union(core_Spleen_Pre, union(core_Spleen_Post1, core_Spleen_Post2)))

# Extract taxonomic information for unique taxa
unique_taxa_post3 <- tax_table(pj2)[unique_post3, ]

# Print unique species
print(unique_taxa_post3)


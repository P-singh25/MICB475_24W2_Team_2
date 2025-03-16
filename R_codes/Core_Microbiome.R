library(phyloseq)
library(microbiome)
library(ggVennDiagram)

sample_data(pj2)

# Subset the phyloseq object for MLN Pre-ICI and MLN Post-ICI3
MLN_Pre <- subset_samples(pj2, location == "MLN" & group == "Pre-ICI")
MLN_Post <- subset_samples(pj2, location == "MLN" & group == "Post-ICI3")

Spleen_Pre <- subset_samples(pj2, location == "Spleen" & group == "Pre-ICI")
Spleen_Post <- subset_samples(pj2, location == "Spleen" & group == "Post-ICI3")

TDLN_Pre <- subset_samples(pj2, location == "TDLN" & group == "Pre-ICI")
TDLN_Post <- subset_samples(pj2, location == "TDLN" & group == "Post-ICI3")

Tumor_Pre <- subset_samples(pj2, location == "Tumor" & group == "Pre-ICI")
Tumor_Post <- subset_samples(pj2, location == "Tumor" & group == "Post-ICI3")


# Transform sample counts to relative abundance (optional but recommended for core analysis)
MLN_Pre_relabund <- transform_sample_counts(MLN_Pre, fun = function(x) x / sum(x))
MLN_Post_relabund <- transform_sample_counts(MLN_Post, fun = function(x) x / sum(x))

Spleen_Pre_relabund <- transform_sample_counts(Spleen_Pre, fun = function(x) x / sum(x))
Spleen_Post_relabund <- transform_sample_counts(Spleen_Post, fun = function(x) x / sum(x))

TDLN_Pre_relabund <- transform_sample_counts(TDLN_Pre, fun = function(x) x / sum(x))
TDLN_Post_relabund <- transform_sample_counts(TDLN_Post, fun = function(x) x / sum(x))

Tumor_Pre_relabund <- transform_sample_counts(Tumor_Pre, fun = function(x) x / sum(x))
Tumor_Post_relabund <- transform_sample_counts(Tumor_Post, fun = function(x) x / sum(x))


# Identify core microbiome for each group
# Core microbiome is defined as taxa present in at least 80% of samples (prevalence = 0.3) with any non-zero abundance (detection = 0)
core_MLN_Pre <- core_members(MLN_Pre_relabund, detection = 0, prevalence = 0.3)
core_MLN_Post <- core_members(MLN_Post_relabund, detection = 0, prevalence = 0.3)

core_Spleen_Pre <- core_members(Spleen_Pre_relabund, detection = 0, prevalence = 0.3)
core_Spleen_Post <- core_members(Spleen_Post_relabund, detection = 0, prevalence = 0.3)

core_TDLN_Pre <- core_members(TDLN_Pre_relabund, detection = 0, prevalence = 0.3)
core_TDLN_Post <- core_members(TDLN_Post_relabund, detection = 0, prevalence = 0.3)

core_Tumor_Pre <- core_members(Tumor_Pre_relabund, detection = 0, prevalence = 0.3)
core_Tumor_Post <- core_members(Tumor_Post_relabund, detection = 0, prevalence = 0.3)


########## Post ICI-3 ############

# Find shared taxa 
shared_taxa_post <- Reduce(intersect, list(core_Tumor_Post, core_MLN_Post, core_TDLN_Post, core_Spleen_Post))
# Prune phyloseq object to shared taxa
shared_physeq_post <- prune_taxa(shared_taxa_post, pj2)
# View taxonomic information
tax_table(shared_physeq_post)

# Plot the venn diagram
core_all_post_ICI3 <- ggVennDiagram(
  x = list(
    "MLN Post-ICI3" = core_MLN_Post,
    "TDLN Post-ICI3" = core_TDLN_Post,
    "Spleen Post-ICI3" = core_Spleen_Post,
    "Tumor Post-ICI3" = core_Tumor_Post
  )
)
core_all_post_ICI3

# Save the venn diagram
ggsave(("Core_microbiome_postICI3.png"), core_all_post_ICI3, width = 10, height = 10)


########## Pre ICI ############

# Find shared taxa
shared_taxa_pre <- Reduce(intersect, list(core_Tumor_Pre, core_MLN_Pre, core_TDLN_Pre, core_Spleen_Pre))
# Prune phyloseq object to shared taxa
shared_physeq_pre <- prune_taxa(shared_taxa_pre, pj2)
# View taxonomic information
tax_table(shared_physeq_pre)


# Plot the venn diagram
core_all_pre_ICI <- ggVennDiagram(
  x = list(
    "MLN Pre-ICI" = core_MLN_Pre,
    "TDLN Pre-ICI" = core_TDLN_Pre,
    "Spleen Pre-ICI" = core_Spleen_Pre,
    "Tumor Pre-ICI" = core_Tumor_Pre
  )
)
core_all_pre_ICI

# Save the venn diagram
ggsave(("Core_microbiome_preICI.png"), core_all_pre_ICI, width = 12, height = 10)



########## Pre ICI and Post ICI3  ############
ggVennDiagram(x=list(core_MLN_Pre, core_MLN_Post))

ggVennDiagram(x=list(core_Spleen_Pre, core_Spleen_Post))

ggVennDiagram(x=list(core_TDLN_Pre, core_TDLN_Post))

ggVennDiagram(x=list(core_Tumor_Pre, core_Tumor_Post))

MICB475: Data Science Research in Microbiology
Team 2 - QingRu Kong, Pranjali Singh, Ran Tao, Tina Wang, Zurui Zhu

Mar 14th, 2025

## Agenda

- Attempted taxonomic bar plots, alpha diversity plots, beta diversity plots, core microbiome analysis, DEseq analysis
- p value should be = or <?

### Taxa Bar Plot: 

Include all locations 
> <img src="../Taxa_Bar.png" height="300">

Spleen only (test)
> <img src="../Spleen_Taxa_Bar.png" height="300">

All locations using another code 
> <img src="../taxonomy_relative_abundance.png" height="300">

### Alpha Diversity 
Using Faith PD
> <img src="../Alpha_Diversity_Pretreatment.png" height="300">

Using Shannon
Based on Organ site
MLN
> <img src="/Figures/Alpha Diversity/Alpha_Diversity_Shannon_MLN.png" height="300">
TDLN
> <img src="/Figures/Alpha Diversity/Alpha_Diversity_Shannon_TDLN.png" height="300">
Spleen
> <img src="/Figures/Alpha Diversity/Alpha_Diversity_Shannon_Spleen.png" height="300">
Tumor
> <img src="/Figures/Alpha Diversity/Alpha_Diversity_Shannon_Tumor.png" height="300">

Based on Treatment timepoint
Pre-ICI
> <img src="/Figures/Alpha Diversity/Alpha_Diversity_Shannon_Pre-ICI.png" height="300">
Post-ICI1
> <img src="/Figures/Alpha Diversity/Alpha_Diversity_Shannon_Post-ICI1.png" height="300">
Post_ICI2
> <img src="/Figures/Alpha Diversity/Alpha_Diversity_Shannon_Post-ICI2.png" height="300">
Post_ICI3
> <img src="/Figures/Alpha Diversity/Alpha_Diversity_Shannon_Post-ICI3.png" height="300">

### Beta Diversity
> <img src="../Beta_Diversity_unweighted_unifrac.png" height="300">

### Core Microbiome 

Pre-ICI 
> <img src="../Pre-ICI for all locations at 0.3 prevelance.png" height="300">

Post-ICI3 
> <img src="../Post-ICI3 for all locations at 0.3 prevelance.png" height="300">

### DEseq: 

all
> <img src="../Deseq_all.png" height="300">

spleen only (Test)
> <img src="../Deseq_Spleen.png" height="300">

#### For next meeting:
- All graphs should be all done by next week


## Meeting Minutes
- Try Kregg for gene analysis
- ok with genus
- detection = 0, prevenenlce = 0.3, approved!
    - evelyn 0.8, ujemi 0.1, ok with in between
 
**Tax bar plot**
- spleen only tax bar plot doesn't have percentage
- should only look at top 1 or the top 1% and other should be other
- axis label too small, (axis.title=) , text size = 16 or 18, title bigger
- theme is not consist with alpha diversity
- need good color for phylum

**Alpha diversity**
- don't need legend for alpha diversity
- use p<
- text ok, title ok, box plot ok, color ok, consistent
- look into false discovery rate, FDR, adjust for increasing false positive, more pairwise = more false positive
- show stat for all three, ns for non significant
- **heat map of significant between conditions, put ANOVA in a heat map**
- Shannon ok, Observed maybe
- use Spleen as "control", compare everything to spleen in the single treatment timepoint

**Beta diversity**
- bit small
- facet differently (2 by 2, not 1 by 4)
- color are not matching alpha diversity, MATCH!!!!!, match theme
- eclipse good
- not really significant for everything even post-ICI3
- include deta diversity, just make an argument for it

**Core Microbiome**
- slightly cut out
- just label organ, don't label treatment condition
  
**DeSeq**
- theme consistent
- change x axis angle at 45 degree
- Do it by organ not all
- 40 is too significant, LOOK INTO THIS!!!
- Tina question: different asv, can't assign to species, just leave as it is, just leave the level;** can blast veill vs veill.1 to compare **
  

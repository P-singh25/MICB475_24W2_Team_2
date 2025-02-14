## MICB475: Data Science Research in Microbiology
Team 2 - QingRu Kong, Pranjali Singh, Ran Tao, Tina Wang, Zurui Zhu

Feb 14th, 2025 

## Agenda

### Research question: 

What is the spatial microbial difference in separate organs at various time points throughout the ICT treatment in mice with melanoma? 

### Aim 1
1. Look at microbiome differences in different organs using Alpha diversity.
2. Weighted UniFrac Diversity analysis (beta diversity) using different grouping methods
   
   a. grouping by organs
   
   b. grouping by timepoints

### Aim 2
1. Track gut bacteria migration through different organs
2. Identify specific population of bacteria that are not previously present in other organs
3. Identify the emergence of new bacteria in the lymphoid organs at the endpoint, and backtrack the origin

### Question for Hans/Evelyn:
1. All the annotations on the figures. Do they have to be using R? Can we annotate using other methods?
2. Dendritic Cell 16 rRNA data was used in the paper, but we could not find it?
3. What can we use to determine metabolic pathways associated with the bacteria we identify? How do we determine if the bacteria is pathogenic?
4. Would we be able to correlate bacteria with disease severity (tumor progression)? (I think we need to find other datasets...)
5. Confirm our truncation length for denoising our data.

### Meeting Minutes
Approach to the paper

Aim 1: compare microbial diversity across organs
1. alpha and beta diversity
2. core microbiome
3. deseq
4. indicator taxa
Aim 2： within each organ site, compare 4 days periods where Day 1-4 is the pre-treatment (control)
1. alpha and beta diversity
2. core microbiome
3. indicator taxa --> will identify microbes of interest
Aim 3: longitudinal of microbes of interest across organs
1. longitudinal analysis
2. optional: network analysis
Aim 4: literature research
1. correlate microbe with ICT/tumor


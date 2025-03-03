# P04 - Aim 1: To determine the alpha diversity across gut, MLN, Spleen, TDLN, and tumor before and after ICI treatment.

Feb 28th, 2025

## Purpose:
To perform the alpha diversity analysis using Faith PD and generate the plot for the alpha diversity analysis. 
- Perform alpha diversity analysis across organs before ICI treatment
- Perform alpha diversity analysis on each organ to compare the changes in diversity throughout the course of ICI treatment

## Material: 
1. R & Rstudio
2. pj2.RData (phyoseq object)

## Method:
1. Import libraries: phyloseq, picante, tidyverse
2. Load the pj2.Rdata<br/>
#### Analysis 1: alpha diversity analysis across organs before ICI treatment
1. Subset the phyloseq object to include only the "Pre-ICI" group -> phylo_pre
2. Perform Faith PD analysis using pd() -> phylo_dist_pre
   - include.root=F
4. Add the Faith PD analysis data (PD) into the phyloseq object (phylo_pre)
5. Plot using ggplot (x = location, y = PD) -> plot.pd_pre
#### Analysis 2: alpha diversity analysis on each organ to compare the changes in diversity throughout the course of ICI treatment
1. Perform Faith PD analysis using pd() on the pj2 phyloseq object
   - include.root=F
3. Add the Faith PD analysis data (PD) into the phyloseq object (pj2)
4. Plot using ggplot (x = location, y = PD) -> plot.pd_location


## Code: 
[Alpha diversity code](../R_codes/Alpha_diversity.R)
   
## Results: 
#### Analysis 1



## Discussion:

## Future direction:


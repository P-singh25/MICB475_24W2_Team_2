# P05 - Aim 2: To identify specific bacterial species that are differentially abundant in the secondary lymphoid organs and tumor post-ICI treatment, and their relationship to the gut microbiome.

March 8st, 2025

## Purpose:
To perform taxa bar plot for each location based on different time groups to determine the changes in relative abundance. 

## Material: 
1. R & Rstudio
2. pj2.RData (phyoseq object)

## Method:
1. Import libraries: phyloseq, ggplot2, dplyr, tidyverse
2. Load the pj2 phyloseq object
3. Melt the pj2 into a dataframe -> pj2_melted, filter the dataframe to remove Day 0 group, and reorder the groups in the order of "Pre-ICI", "Post-ICI1", "Post-ICI2", "Post-ICI3", remove NA.
4. Calculate the total number of ASVs for each group of each location, generate the new column Total_Abundance
5. Calculate the relative abundance for each group of each location, generate the new column Relative_Abundance
6. Use ggplot to plot the relative abundance plot (x = group, y = Relative_Abundance, fill = phylum, facet_wrap based on location)

## Code: 
[Taxonomic Bar Plot code](../R_codes/"Taxonomic_Bar_Plot_FINAL.R") 
   
## Results: 
#### 



## Discussion:

## Future direction:

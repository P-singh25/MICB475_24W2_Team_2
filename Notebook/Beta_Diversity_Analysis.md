# P06- To determine changes in microbial composition over ICI treatment course in spleen, MLN, TDLN and tumor

March 8st, 2025

## Purpose:
To perform beta diversity (Bray Curtis) analysis to plot the abundance of each sample to determine changes in microbial composition across location and over time.

## Material: 
1. R & Rstudio
2. pj2.RData (phyoseq object)

## Method:
1. Import libraries: tidyverse, phyloseq, vegan, and ggplot2
2. Filter phyloseq object to keep group: Pre-ICI, Post-ICI1, Post-ICI2, Post-ICI3, and location: spleen, MLN, TDLN, tumor.
3. Load phyloseq object after filtering and rarification
4. Convert ASV matrix from phyloseq object into a dataframe
5. Convert metadata from phyloseq object into a dataframe
6. Perform Bray Curtis Principal Coordinate Analysis and plot using ggplot2. Draw eclipses to represent 95% confidence interval of the each condition. PERMANOVA statistical analysis test to generate R-squared statistic and p-value in a results table.

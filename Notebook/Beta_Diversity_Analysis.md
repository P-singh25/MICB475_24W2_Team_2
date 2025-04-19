Method:
Import libraries: tidyverse, phyloseq, vegan, and ggplot2
Filter phyloseq object to keep treatment type == Separate and experiment group == general and cohouse
Load phyloseq object after filtering and rarification
Convert ASV matrix from phyloseq object into a dataframe
Convert metadata from phyloseq object into a dataframe
Filter the metadata to only keep variables that are related to gut microbiome dysbiosis(Sample.Name, Age.New.Bin, Cage.ID, Experiment.Group, Genotype, Mouse.ID, Phenotype.score, FD.severity, Sex, Weight.grams)
Create For Loops to iterate over each variable in the metadata, remove missing data, filter out NA values, calculate a dissimilarity matrix based on Bray-Curtis distance, and perform the PERMANOVA statistical analysis test, controlling for cage ID, to generate R-squared statistic and p-value in a results table.
Adjust p-values
Filter the result table to include only significant variables with Padjust < 0.05
Generate a bar plot using ggplot2 to visualize the R-squared values for each significant variable

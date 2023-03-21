Authors: Shivkumar Vishnempet Shridhar†, Francesco Beghini†, Marcus Alexander, Ilana L. Brito*, and Nicholas A. Christakis*
By: Christakis Group (Yale), and Brito Group (Cornell), United States of America

## Honduras Microbiome phenotype project

This github repo describes workflow and codes used in Honduras microbiome study:

### Contents:

- Species and pathway abundance profiles
- Microbiome-phenotype association 
- Microbiome-pathway association 
- Calculation of microbiome variance explained by phenotypes and pathways
- Plotting scripts

### Species abundance profiles

Relative abundances of species and pathways generated by Metaphlan4 and Humann3 respectively after pre-processing.

### Microbiome-phenotype association

We evaluated associations between gut microbiome and phenotypes using linear mixed model.

  Species abundance ~ age + sex + BMI + batch effect + bristol stool scale + DNA concentration + Sampling date + 1|village + phenotype
  
where the relative abundance of species was transformed using the centred additive log-ratio (CLR) transformation.

### Microbiome-pathway association

We evaluated associations between metabolic pathways and phenotypes using linear mixed model.

  Pathway abundance ~ age + sex + BMI + batch effect + bristol stool scale + DNA concentration + Sampling date + 1|village + phenotype

### Calculation of microbiome variance explained by phenotypes and pathways

The microbiome composition variance explained by phenotypes was calculated by permutational multivariate analysis of variance using distance matrices, implemented in the adonis function for R package vegan (v.2.6), using 1000 permutations and a Bray-Curtis distance matrix calculated using relative abundances of microbial species. Variance explained was also performed using relative abundances of MetaCyc microbial biochemical pathways separately.

### Plotting scripts

Consisting of plotting scripts for all figures


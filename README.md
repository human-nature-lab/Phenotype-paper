Authors: Shivkumar Vishnempet Shridhar†, Francesco Beghini†, Marcus Alexander, Ilana L. Brito*, and Nicholas A. Christakis*
By: Christakis Group (Yale), and Brito Group (Cornell), United States of America

## Honduras Microbiome phenotype project

This github repo describes workflow and codes used in Honduras microbiome study:

### Contents:

- Species abundance profiles
- Microbiome-phenotype association 
- Microbiome-pathway association 
- Calculation of microbiome variance explained by phenotypes and pathways
- Plotting scripts

### Species abundance profiles

Metagenomes were profiled consistent with previous data analysis of 1000IBD and Lifelines-DEEP10 cohorts, as follows. KneadData tools (v0.5.1) were used to process metagenomic reads (in fastq format) by trimming the reads to PHRED quality 30 and removing Illumina adapters. Following trimming, the KneadData integrated Bowtie2 tool (v2.3.4.1) was used to remove reads that aligned to the human genome (GRCh37/hg19).
Taxonomic composition of metagenomes was profiled by MetaPhlAn2 tool (v2.7.2) using the MetaPhlAn database of marker genes mpa_v20_m200. Profiling of genes encoding microbial biochemical pathways was performed using the HUMAnN2 pipeline (v0.11.1) integrated with the DIAMOND alignment tool (v0.8.22), UniRef90 protein database (v0.1.1) and ChocoPhlAn pan-genome database (v0.1.1). As a final quality control step, samples with unrealistic microbiome composition (eukaryotic or viral abundance > 25% of total microbiome content or total read depth < 10 million) were excluded, leaving 8,208 samples for further analyses. Analyses were performed using locally installed tools and databases on CentOS (release 6.9) on the high-performance computing infrastructure available at our institution and using the MOLGENIS data platform2.

Rarefaction and extrapolation (R/E) sampling curves for estimation of total richness of species and genera in the population were constructed using a sample size-based interpolation/ extrapolation algorithm implemented in the iNEXT package for R. Variance of phyla in the population was compared to variance of pathways using F Test for comparison of variances and total variance of microbial taxa vs pathways were compared using permutation tests. 

- An example of metagenome processing jobs is provided in *microbiome_profiling* folder
- Rarefaction curves, variance comparison and other microbiome description codes are in *microbiome_description* folder

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


# Transcriptonal Analysis of Stress Responses
The files here were part of my research under Dr. William Kath at the NSF-Simons Center for Quantitative Biology. The research term was 3 months long, ending in August 2021.

The goal of the project was to investigate multi-kinase control of environmental stress responsive transcription in yeast. What kinases transmit stress signals and what roles of theirs are contingent on the environment? My research sought to build on research presented in Mace PLOS 2020 that showed how kinases are involved in regulating genes in particular environments.


## DESeq Data Processing and Exploration

I aligned DESeq transcripts to a yeast reference genome using STAR in a computing cluster. The data and alignment scripts are not included here. The stresses on each run are listed in Pincus_SraRunTable.csv in addition to information about each run's varied perturbation and mutated genotype. The file count_tables.R creates a count table from the transcripts in a DESeq matrix. I familiarize myself with the counts data by visualizing the DESeq data with heatmaps and looking at the most differently expressed genes in the dataset. 


## Dimensionality Reduction

The files yeast_all_PCA_v2.R and multik_counts_all.R reduce the transcripts via Principal Component Analysis. Number of PCs was selected with bootstrapping performed in multik_boot.R and shown in results/boot_plot.pdf. Further analysis used the first 25 components that accounted for most of the variation in the counts data. PCA loadings can be viewed by running loadings.R.

Figures showing runs plotted in PC1-PC2 space and colored by genotype or stress type can be found in results/genotype_scaled.pdf and results/multik_stress_types.pdf. The horizontal clustered structure revealed in the PCA plots as well as the t-SNE plot in results/tsne.pdf suggest a patterned response to environmental stress that was further investigated later on. 

Using GO analysis and further research on the genotypes clustered together in PCA-reduced space, I reasoned that the structure revealed in the PCA plot of the first 2 components forms due to common nucleic acid binding activity. 



## Resources
### Mace PLOS 2020
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0230246




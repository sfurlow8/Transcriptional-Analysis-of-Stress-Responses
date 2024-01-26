setwd("/Users/sfurlow/Desktop/CQuB/alignment_pca")


library(ggplot2)
library(ggcorrplot)
library(GGally)

gene_counts <- read.table(file="yeast_WT_SNF2_all_counts.tsv",sep="",header=TRUE)
head(gene_counts)

rownames(gene_counts)<-gene_counts$gene
gene_counts$gene<-NULL

sample_names<-names(gene_counts)
sample_types<-matrix(unlist(strsplit(sample_names,"_")),ncol=2,byrow=TRUE)[,1]

# experiment table links each sample to their condition (WT or SNF)
expt_table<-data.frame(sample_names=sample_names,sample_types=as.factor(sample_types))

library(DESeq2)

# DESeq standardizes counts
dds <- DESeqDataSetFromMatrix(countData=gene_counts,colData=expt_table,design=~sample_types)
#rld <- rlog(dds, blind=FALSE)
vsd <- vst(dds, blind=FALSE) # removes dependence of variance on mean: 
# otherwise highly exp. genes would dominate determination of PCA directions

vsd2<-assay(vsd)[which(apply(assay(vsd),1,FUN=var)>0),]

# samples as rows, genes as columns
# compare correlations between genes
gene_pca<-prcomp(t(vsd2), center=TRUE, scale=FALSE)

summary(gene_pca) # look at principal components
screeplot(gene_pca,type="l", npcs=15) # use elbow rule: 5 or 6 PCs

percentVar <- gene_pca$sdev^2 / sum( gene_pca$sdev^2 ) # pct of total var in each component

intgroup.df <- as.data.frame(colData(vsd)[,"sample_types", drop=FALSE])
group<-colData(vsd)[["sample_types"]]

# df with PC1, PC2, condition, and name of each sample
d <- data.frame(PC1=gene_pca$x[,1], PC2=gene_pca$x[,2], group=group, intgroup.df, name=colnames(vsd))

ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
  coord_fixed(ratio=sqrt(percentVar[2]/percentVar[1])) +
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) 

library(PCAtools)
p <- pca(vsd2, metadata = colData(dds), removeVar = 0.05)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p, showLoadings = TRUE, labSize = 5, pointSize = 5, sizeLoadingsNames = 2, ntopLoadings = 10, maxoverlapsConnectors = 25)
plotloadings(p, labSize = 2)

gene_loadings<-p$loadings
pca_load1<-gene_loadings[,1,drop=FALSE]
o <- order(abs(pca_load1$PC1), pca_load1$PC1, decreasing = TRUE)
pca_load1_sorted<-pca_load1[o,,drop=FALSE]
pca_load2<-gene_loadings[,2,drop=FALSE]
o <- order(abs(pca_load2$PC2), pca_load2$PC2, decreasing = TRUE)
pca_load2_sorted<-pca_load2[o,,drop=FALSE]

head(pca_load1_sorted,20)
head(pca_load2_sorted,20)

pca_genes<-gene_pca$rotation
pca_gene1<-pca_genes[,1,drop=FALSE]
o <- order(abs(pca_gene1), decreasing = TRUE)
pca_gene1_sorted<-pca_gene1[o,1,drop=FALSE]
head(pca_gene1_sorted,20)
pca_gene2<-pca_genes[,2,drop=FALSE]
o <- order(abs(pca_gene2), decreasing = TRUE)
pca_gene2_sorted<-pca_gene2[o,1,drop=FALSE]
head(pca_gene2_sorted,20)


library(ggplot2)
library(ggcorrplot)
library(GGally)

setwd("/Users/sfurlow/Desktop/CQuB/alignment_pca")


## CREATE A COUNT TABLE GENE_COUNTS
gene_counts <- as.data.frame(matrix(nrow=7131, ncol=12))

for (ii in 1:6) {
  str_wt = paste('yeast_0', ii, '_WT_', 'ReadsPerGene.out.tab', sep='')
  str_snf = paste('yeast_0', ii, '_SNF2_', 'ReadsPerGene.out.tab', sep='')
  colnames(gene_counts)[ii] = paste('WT_0', ii, sep='') # name columns
  colnames(gene_counts)[ii+6] = paste('SNF_0', ii, sep='')
  
  wt_data  <- read.table(str_wt, sep="", header=FALSE)
  gene_counts[,ii] <- wt_data[,4] # 3rd column for pair ended reads
  snf_data <- read.table(str_snf, sep="", header=FALSE)
  gene_counts[,ii+6] <- snf_data[,4] 
}
gene_counts <- gene_counts[-c(1:4),] # take genes only
rownames(gene_counts) = wt_data[-c(1:4),1]

gene_counts <- gene_counts[,-c(6,12)] # remove sample 6 WT and SNF
gene_counts <- t(gene_counts)

## WORKING WITH DESEQ
samplename<-colnames(gene_counts)
samplematrix<-matrix(unlist(strsplit(samplename,'_')),ncol=2,byrow=TRUE)
sample_condition<-samplematrix[,1]
exp_matrix<-data.frame("sample"=samplename,"condition"=sample_condition)
DDS_yeast <- DESeqDataSetFromMatrix(countData=gene_counts,colData=exp_matrix,design=~condition)


DDS_yeast_result <- DESeq(DDS_yeast)
yeast_normed <- counts(DDS_yeast_result,normalized=TRUE)
write.csv(as.data.frame(yeast_normed),file="yeast_normed.csv")

yeast_results <- results(DDS_yeast_result, lfcThreshold=0.5, alpha=0.05,
                         altHypothesis="greaterAbs", contrast=c("condition","SNF","WT"))
yeast_ordered <- yeast_results[order(yeast_results$padj),]
write.csv(as.data.frame(yeast_ordered),file="yeast_ordered.csv")

# yeast_combined.csv contains all of the results in one file, 
# with genes ordered by the adjusted p-value associated with the log2 fold change
yeast_merged<-merge(yeast_normed,yeast_ordered,by="row.names",sort=FALSE)
rownames(yeast_merged)<-yeast_merged$Row.names
yeast_merged$Row.names<-NULL
yeast_combined <- yeast_merged[order(yeast_merged$padj),]
write.csv(as.data.frame(yeast_combined),file="yeast_combined.csv")

# visualize DESeq results
plotMA(yeast_results, ylim=c(-2.5,2.5)) # plotMA summarizes fold change results
yeastLFC <-lfcShrink(DDS_yeast_result, coef="condition_WT_vs_SNF", lfcThreshold =0.5)
plotMA(yeastLFC , ylim=c( -2.5 ,2.5))


## EXPLORING THE DATA WITH HEATMAPS
vsd  <- vst(DDS_yeast_result , blind=FALSE)
library("pheatmap")

# sorts the rows in the count table by decreasing mean expression and then pulls out the top 20
select  <- order(rowMeans(counts(DDS_yeast_result ,normalized=TRUE)),
                 decreasing=TRUE)[1:20] 

df <- as.data.frame(colData(DDS_yeast_result)[,c("condition","sample")])

# heatmap gives some idea of how the samples are clustered
pheatmap(assay(vsd)[select ,],  cluster_rows=TRUE , show_rownames=TRUE ,
         cluster_cols=TRUE , annotation_col=df)



## PRINCIPAL COMPONENT ANALYSIS
library(ggplot2)

pcaData  <- plotPCA(vsd , intgroup=c("condition", "sample"), returnData=TRUE)
percentVar  <- round (100 * attr(pcaData , "percentVar"))

ggplot(pcaData , aes(PC1 , PC2 , color=sample , shape=condition)) +
  coord_fixed(ratio = sqrt(percentVar [2]/ percentVar [1])) +
  geom_point(size =3) +
  xlab(paste0 ("PC1: ",percentVar [1],"%  variance")) +
  ylab(paste0 ("PC2: ",percentVar [2],"%  variance"))


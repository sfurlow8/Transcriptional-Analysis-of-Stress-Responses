setwd("/Users/sfurlow/Desktop/CQuB/alignment_pca")

library(ggplot2)
library(ggcorrplot)
library(GGally)



## PROCESSING-----------------------

# load count table
gene_counts = read.table("/Users/sfurlow/Desktop/count_table.txt", header=TRUE)
head(gene_counts)

# rename and reorder the samples columns in gene_counts
for (ii in 1:length(colnames(gene_counts))) {
  if (startsWith(colnames(gene_counts)[ii], "X00")) {
    suff = unlist(strsplit(colnames(gene_counts)[ii], "X"))[2]
    colnames(gene_counts)[ii] = paste("S8", suff, sep='')
  }
  else {
    suff = unlist(strsplit(colnames(gene_counts)[ii], "X"))[2]
    colnames(gene_counts)[ii] = paste("S7", suff, sep='')
  }
}
gene_counts=gene_counts[,order(colnames(gene_counts))]
head(gene_counts)

# experiment matrix
expt_mat = read.csv("Pincus_SRARunTable.csv", header=TRUE)
rownames(expt_mat)<-expt_mat$Run
expt_mat$Run<-NULL
# expt_mat[,-c(1)]=apply(expt_mat[,-c(1)], MARGIN = 2, FUN = factor)
expt_mat$Genotype = factor(expt_mat$Genotype) # factor the conditions
expt_mat$kinase_inhibitor = factor(expt_mat$kinase_inhibitor)
expt_mat$media = factor(expt_mat$media)
expt_mat$perturbation = factor(expt_mat$perturbation)
expt_mat$stress = factor(expt_mat$stress)
head(expt_mat)

expt_mat = expt_mat[order(rownames(expt_mat)),] # reorder rows to match order of gene_counts
head(expt_mat)

# correlation plot
C<-cor(gene_counts)
p<-ggcorrplot(C,tl.cex=5)
p+scale_fill_gradient2(limit=c(0.9,1),low="blue",high="red",mid=" white",midpoint =0.950)


## DESEQ TRANSFORMATIONS-------------------

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=gene_counts,colData=expt_mat,design=~1)

vsd <- vst(dds, blind=FALSE) # variance stabilizing transformation
vsd2<-assay(vsd)[which(apply(assay(vsd),1,FUN=var)>0.05),]



## PCA & SCREE PLOTS--------------------------

# samples as rows, genes as columns
# compare correlations between genes
gene_pca<-prcomp(t(vsd2), center=TRUE, scale=FALSE)
gene_pca<-prcomp(t(vsd2), center=TRUE, scale=TRUE) # scaling

summary(gene_pca) # look at principal components
screeplot(gene_pca,type="l", npcs=15) # uninstall PCAtools to run this

#library(factoextra)
#fviz_eig(gene_pca, main="Genotype Scree Plot") # screeplot barplot

library(PCAtools)

p <- pca(vsd2, metadata = colData(dds), removeVar = 0.05)

elbow=findElbowPoint(p$variance) # elbow heuristic
min=which(cumsum(p$variance) > 75)[1] # min components to cover certain cumulative var

screeplot(p, components=getComponents(p)[1:15], title = 'SCREE Plot', 
          vline = c(min, elbow)) +
          geom_label(aes(x = min + 1, y = 50,
                      label = 'CumSum', vjust = -1, size = 8)) +
          geom_label(aes(x = elbow + 1, y = 50,
                      label = 'Elbow', vjust = -1, size = 8))



## BIPLOTS & PCA PLOTS---------------------------

biplot(p, showLoadings = TRUE, labSize = 5, pointSize = 5, 
       sizeLoadingsNames = 2, ntopLoadings = 10, maxoverlapsConnectors = 25)
plotloadings(p, labSize = 2)

percentVar <- gene_pca$sdev^2 / sum( gene_pca$sdev^2 ) # pct of total var in each component
intgroup.df <- as.data.frame(colData(vsd)[,"Genotype", drop=FALSE]) # gives each sample and correspoding group label
group<-rep("Other", 391)
group[colData(vsd)[["Genotype"]]=="HOG1-as"]<- "HOG1-as"
group[colData(vsd)[["Genotype"]]=="PBS2-as"]<- "PBS2-as"
group[colData(vsd)[["Genotype"]]=="STE11-as"]<- "STE11-as"
group[colData(vsd)[["Genotype"]]=="TPK123-as"]<- "TPK-123"


group<-colData(vsd)[["Genotype"]] # gives levels of the groups (group is a factor)
group<-colData(vsd)[["stress"]] # gives levels of the groups (group is a factor)


# PC1 and PC2 in d give the coordinates of each sample for the first 2 components
# d also gives group level and sample name
d <- data.frame(PC1=gene_pca$x[,1], PC2=gene_pca$x[,2], PC3=gene_pca$x[,3], PC4=gene_pca$x[,4],
                PC5=gene_pca$x[,5], PC6=gene_pca$x[,6], group=group, intgroup.df, name=colnames(vsd))

ggplot(data= d %>% filter(group=="Other"), aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=2)+
  geom_point(data = d %>% filter(group!="Other"), size=3) +
  scale_color_manual(values = c("blue", "black", "red", "green", "yellow")) +
  coord_fixed(ratio=sqrt(percentVar[2]/percentVar[1])) +
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  ggtitle("Genotype groups with scaling")


ggplot(data= d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) +
  coord_fixed(ratio=sqrt(percentVar[2]/percentVar[1])) +
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  ggtitle("Stress groups with scaling")


# 3D plot with plotly
library(plotly)
ar21<-sqrt(percentVar[2]/percentVar[1])
ar31<-sqrt(percentVar[3]/percentVar[1])
fig <- plot_ly(d, x = ~PC1, y = ~PC2, z = ~PC3, color = ~group)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC2'),
                                   zaxis = list(title = 'PC3'),
                                   aspectmode='manual',
                                   aspectratio = list(x=1, y=ar21, z=ar31)))

fig


tester_count = gene_counts
none_kin_count = tester_count[,which(tester$kinase_inhibitor=="None")]
tester_exp = expt_mat
none_kinase = tester[which(tester$kinase_inhibitor=="None"),]


# t-SNE plot
library(M3C)

tsne(t(gene_counts),labels=as.factor(expt_mat$stress), perplex=100,
     legendtitle = "Stress")


library(Rtsne)
set.seed(8)  
tsne_model_1 = Rtsne(t(as.matrix(vsd2)), check_duplicates=TRUE, pca=TRUE, perplexity=100, theta=0.0, dims=2, max_iter = 5000)

# data frame from sra run table
sra_file<-read.table(file="/Users/sfurlow/Downloads/SraRunTable1.txt",sep=",",header=TRUE)
uniques<-apply(sra_file,2,unique)

sample_data<-data.frame(cbind(sra_file$Run,sra_file$Genotype,sra_file$kinase_inhibitor,sra_file$media,sra_file$perturbation,sra_file$stress))
colnames(sample_data)<-c("Run","Genotype","kinase_inhibitor","media","perturbation","stress")
rownames(sample_data)<-paste0("X",substr(sample_data$Run,8,10))

## getting the two dimension matrix
d_tsne_1 = as.data.frame(tsne_model_1$Y)
sample_names<-rownames(expt_mat)
d_tsne_1$Genotype<-as.factor(sample_data$Genotype)
d_tsne_1$kinase_inhibitor<-as.factor(sample_data$kinase_inhibitor)
d_tsne_1$media<-as.factor(sample_data$media)
d_tsne_1$perturbation<-as.factor(sample_data$perturbation)
d_tsne_1$stress<-as.factor(sample_data$stress)
#rownames(d_tsne_1)<-rownames(TPMcoldata)
#d_tsne_1$type<-TPMcoldata$type
#d_tsne_1$sex<-TPMcoldata$sex
#d_tsne_1$sample<-TPMcoldata$sample

ggplot(d_tsne_1, aes(x=V1, y=V2, colour=stress)) +
  geom_point(size=3)

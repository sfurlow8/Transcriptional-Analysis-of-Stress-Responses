pca_genes<-gene_pca$rotation
pca_gene1<-pca_genes[,1,drop=FALSE] # loadings for PC1
o <- order(abs(pca_gene1), decreasing = TRUE) # indices of ordered loadings
pca_gene1_sorted<-pca_gene1[o,1,drop=FALSE] # ordered loadings
cat(head(rownames(pca_gene1_sorted),50))

pca_gene2<-pca_genes[,2,drop=FALSE]
o <- order(abs(pca_gene2), decreasing = TRUE)
pca_gene2_sorted<-pca_gene2[o,1,drop=FALSE]
cat(head(rownames(pca_gene2_sorted),50))

pca_gene5<-pca_genes[,5,drop=FALSE]
o <- order(abs(pca_gene5), decreasing = TRUE)
pca_gene5_sorted<-pca_gene5[o,1,drop=FALSE]
cat(head(rownames(pca_gene5_sorted),50))
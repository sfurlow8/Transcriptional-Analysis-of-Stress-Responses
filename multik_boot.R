setwd("/Users/sfurlow/Desktop/CQuB/alignment_pca")

## BOOTSTRAP 

boot_eigs = as.data.frame(matrix(nrow=ncol(gene_counts), ncol=100))
rownames(boot_eigs) = paste0("eig",rownames(boot_eigs))

for (ii in 1:50) {
  samps = apply(gene_counts, MARGIN=c(1,2), FUN=sample,
                size=nrow(boot_eigs), replace=TRUE)
  samps_pca = prcomp(samps, scale. = TRUE)
  boot_eigs[,ii] = (samps_pca$sdev)^2
  
}
avg_boot_eigs = apply(boot_eigs, MARGIN=1, FUN=mean)
head(avg_boot_eigs)


##---------------------------------------------------


nits=100
nsamp=ncol(vsd2)
evsum<-rep(0,nsamp)
for (k in 1:nits) {
  print(k)
  newcols<-sample(1:nsamp,replace=TRUE) # new column sample indices
  vsd_tmp<-vsd2[,newcols] # new sample w replacement
  pca_tmp<-prcomp(t(vsd_tmp),center=TRUE,scale=TRUE) # pca on t(vsd2)
  ev_tmp<-pca_tmp$sdev^2 
  evsum<-evsum+ev_tmp # total eig sum
}
evsum<-evsum/nits # avg eigs

evnull<-rep(0,nsamp)
for (k in 1:nits) {
  print(k)
  vsd_tmp<-t(apply(vsd2,1,sample,replace=TRUE))
  pca_tmp<-prcomp(t(vsd_tmp),center=TRUE,scale=TRUE)
  ev_tmp<-pca_tmp$sdev^2
  evnull<-evnull+ev_tmp
}
evnull<-evnull/nits # avg eigs

plot(evsum,xlim=c(0,100),ylim=c(1,50))
#lines(ev)
lines(evnull,col="red")

pcs=rep(1:391, 2)
type=c(rep("boot", 391), rep("null", 391))
dat = data.frame(PC=pcs, eigs=c(evsum,evnull), type=type)

ggplot(data=dat, aes(x=PC, y=eigs, colour=type)) + 
  geom_point(data=dat[dat$type=="boot", ],size=3) + 
  geom_line(size=2) +
  xlim(c(0,100)) +
  ylim(c(0,100)) +
  xlab("Number of PCs") +
  ylab("Eigenvalue") +
  theme(text = element_text(size=20))
#,
    #    axis.text.x = element_text(angle=90, hjust=1)) 

ggplot(data=dat, aes(x=PC, y=eigs, color=type)) + geom_point(size=2) + geom_line()+ 
  xlim(c(0,50)) + 
  ylim(c(0, 35)) + # include ylim when you want to zoom in
  xlab("Number of PCs") + ylab("Eigenvalue") + 
  ggtitle("Bootstrapped Eigenvalues")


##---------------------------------------------------


## NULL DATA SET WITH MARGINAL RESAMPLING

null_eigs = as.data.frame(matrix(nrow=ncol(gene_counts), ncol=100))
rownames(null_eigs) = paste0("eig",rownames(null_eigs))

# don't use raw counts!! 
# biggest variation comes from differences between read sizes, a few will dominate
# would have to try DESeq normalized counts if you want to use counts
# genes that are highly expressed dominate -> log transform

# generate 100 matrices 
for (jj in 1:100) {
  samps = apply(gene_counts, MARGIN=1, FUN=sample, 
        size=ncol(gene_counts), replace=TRUE)
  mat = t(samps) # mat has genes as rows, 391 columns
  mat_pca = prcomp(mat, scale. = TRUE)
  null_eigs[,jj] = (mat_pca$sdev)^2 # put 391 eigenvalues into null_eigs
}

# across 100 shuffled samples (columns of mat), find the average of each eigenvalue
# avg val for eig1, eig2, eig3, etc
avg_null_eigs = apply(null_eigs, MARGIN=1, FUN=mean)
head(avg_null_eigs)


# experiment with boot function
# have to make a statistic function to input to boot
# pca_eigs function takes a df and indices, returns eigenvalues of pca performed on subset of df
library(boot)
pca_eigs<-function(df, freq) {
  df_sample = df[freq,]
  pca_df = prcomp(df_sample, scale. = TRUE)
  eigs = (pca_df$sdev)^2
  return(eigs)
}

boot_dat=boot(data=gene_counts, statistic=pca_eigs, R=100, sim="parametric")

# plot boot and null eigenvalues vs. # PCs
pcs=rep(1:391, 2)
type=c(rep("boot", 391), rep("null", 391))
dat = data.frame(PC=pcs, eigs=c(boot_dat$t0,avg_null_eigs), type=type)

ggplot(data=dat, aes(x=PC, y=eigs, color=type)) + geom_point(size=2) + geom_line()+ 
  xlim(c(0,50)) + 
  ylim(c(0, 35)) + # include ylim when you want to zoom in
  xlab("Number of PCs") + ylab("Eigenvalue") + 
  ggtitle("Bootstrapped Eigenvalues")



#!/usr/bin/env Rscript

## Aim : we try different window size for DESeq2 analysis. Here, we make a  histogram for 100 and 2048 window size analyses.
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+


######################################################
##
## Look at results from DESeq with different cutoff
##
######################################################

siteSize=2048
treatment='Copper'
strand='both'
window.size=100
numSam = 6
all.name = paste0(treatment,".", siteSize, ".", strand)


out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".alt.run/output/DESeq2.pval.discrete.",10,".Robj")

load(out.path)




#################################################
## make a histogram of statistics from null and alternative and p-values for DESeq
#################################################

setwd("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/fig/DESeqdebug/")
pdf("hist.statistic.pval.DESeq2.100.10.discrete.pdf")

#### filter.cut = 0 window 100

filter.cut = 10
pval = pval.deseq.100.10

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".alt.run/output/DESeq2.min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".null.run/output/DESeq2.min.pval.", filter.cut, ".txt"))[,1])

 
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))

xmax = 1
xmin = 0

par(mfrow = c(3,1))
hist(deseq.alt[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : alt test statistic [100]"), xlim=c(xmin, xmax), xlab ="p-value")
hist(deseq.null[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : null test statistic [100]"), xlim=c(xmin, xmax), xlab ="p-value")
hist(pval[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : pvalue [100]", length(pval)-length(del.ix.deseq)), xlim=c(xmin, xmax), xlab ="p=value")



dev.off()








#################################################
## make a histogram of statistics from null and alternative and p-values for DESeq (pooling)
#################################################

setwd("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/fig/DESeqdebug/")
pdf("hist.statistic.pval.DESeq2.100.10.discrete.pooling.pdf")



#### filter.cut = 0

filter.cut = 10


deseq.alt = as.numeric(as.matrix(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".alt.run/output/DESeq2.all.pval.", filter.cut, ".txt"))))
deseq.null = as.numeric(as.matrix(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".null.run/output/DESeq2.all.pval.", filter.cut, ".txt"))))

del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))


xmax = 1
xmin = 0

par(mfrow = c(2,1))
hist(deseq.alt[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : alt test statistic [100]"), xlim=c(xmin, xmax), xlab ="p-value")
hist(deseq.null[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : null test statistic [100]"), xlim=c(xmin, xmax), xlab ="p-value")


dev.off()







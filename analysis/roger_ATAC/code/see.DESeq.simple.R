#!/usr/bin/env Rscript

## Aim : Resutls from DESeq analysis are very bad. I'm not sure I used DESseq properly. So I'll try different ways to use DESeq and see if results make sens. For this purpose, I'll focus on `Copper.2048.both.300'. I'll try to exclude 300bp windows with low read counts.
## Here, I particularly focus on one filtering (at least 60 reads) and compared results from DESeq and DESeq2. 
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+


wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/'
siteSize=2048
treatment='Copper'
null=TRUE
strand='both'
#strand='plus'
#strand='minus'
window.size=300
numSam = 6


library("DESeq")
library("DESeq2")


## set up working directory 
setwd(wd.path)

## make directory name  
if(!null){
    dir.name=paste0(treatment,".", siteSize, ".", strand,  ".", window.size, ".alt")
}else{
    dir.name=paste0(treatment,".", siteSize, ".", strand, ".", window.size, ".null")
}

## get number of sites informaiton 
path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/", treatment, ".", siteSize, ".numSites.txt")
numSites.list = scan(path)

## read data and check if data is proper
input.dir.path = paste0(wd.path, dir.name, ".run/")
numC = siteSize%/%window.size
deseq.data = matrix(data=NA, nc = numSam, nr = sum(numSites.list)*numC)
st.chr.ix = NULL
en.chr.ix = NULL
st.chr.ix[1] = 1
err.check = rep(NA, 22)
for(chr in 1:22){
    #chr = 1
    numSites = numSites.list[chr]
    numRow = numSites*numC
    en.chr.ix[chr] = st.chr.ix[chr] + numRow - 1
    deseq.data.each = read.table(paste0(input.dir.path, "data.", chr, ".txt"))
    if(dim(deseq.data.each)[1] == numRow){
        err.check[chr] = FALSE
        deseq.data[st.chr.ix[chr]:en.chr.ix[chr],] = as.matrix(deseq.data.each)
    }
    st.chr.ix[chr+1] = en.chr.ix[chr] + 1
}
#sum(err.check==FALSE)
# 22
#sum(is.na(deseq.data))
# 0

#########################
## For DESeq analyses
#########################

## prepare input format for deseq
condition=factor(c(rep("treated", numSam/2),rep("untreated", numSam/2)))
## prepare data 
deseq.full.data =newCountDataSet(deseq.data,condition)


## let's try 60, 30, 20, 10

rsum = rowSums ( counts ( deseq.full.data ))
filter.cut = 60
##use = ((rsum > filter.cut) & (rsum < 6500))
use = (rsum > filter.cut)

deseq.data.filtered = deseq.full.data[ use, ]

## estimate size factors
deseq.data.filtered = estimateSizeFactors(deseq.data.filtered)
sizeFactors(deseq.data.filtered)

## estimate dispersion parameters
deseq.data.filtered=estimateDispersions(deseq.data.filtered)

## perform test
resDESeq.filtered =nbinomTest(deseq.data.filtered,"treated", "untreated")
pval.vec = rep(NA, length(use))
pval.vec[use==TRUE] = resDESeq.filtered$pval
pval.filtered =matrix(pval.vec,ncol=numC,byrow=T)


DESeq.simple.pval = resDESeq.filtered$pval


##########################
## For DESeq2 Analyses
##########################


countData = counts ( deseq.full.data )

###Converts the colData table to a dataframe with vector labels
colData <- data.frame(row.names= c("T1","T2", "T3", "C1", "C2", "C3") ,t=c("T","T","T","C","C","C"),r=c(1,2,3,1,2,3))
### should be a factor!


rsum = rowSums (countData)
filter.cut = 60
use = (rsum > filter.cut)

countData.filtered = countData[ use, ]

###Perform DESeq2 Analysis

ddsTvC <- DESeqDataSetFromMatrix(
	countData=countData.filtered,
	colData=colData,
	design = ~ r + t)
dds <- DESeq(ddsTvC)
res <- results(dds)

DESeq2.simple.pval = res$pvalue


##############################
# QQ plot and Histgram
##############################

setwd("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/fig/DESeqdebug/")

png(paste("DESeq",".qqPlot.png",sep=""))
qqplot(-log10(ppoints(length(DESeq.simple.pval))),-log10(DESeq.simple.pval),xlab="Expected -log10(Pvalue)",ylab="Observed -log10(Pvalue)")
abline(0,1,col='red')
dev.off()

png(paste("DESeq2",".qqPlot.png",sep=""))
qqplot(-log10(ppoints(length(DESeq2.simple.pval))),-log10(DESeq2.simple.pval),xlab="Expected -log10(Pvalue)",ylab="Observed -log10(Pvalue)")
abline(0,1,col='red')
dev.off()

png(paste("DESeq2andDESeq",".Hist.png",sep=""))
par(mfrow =c(2,1))
hist(DESeq.simple.pval, breaks=1000)
hist(DESeq2.simple.pval, breaks=1000)
dev.off()



#############################################################
# Conclusion: My DESeq2 analyses results are similar to what Roger has!
#############################################################









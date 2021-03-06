#!/usr/bin/env Rscript

## Aim : compute empirical p-values for DESeq2 100bp window and full window analyses
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+




##################################################################
# get p-value from empirical null distribution of test statistic (DESeq2)
##################################################################



get.pval.from.empirical.null.dist.discrete <- function(statistic.null, statistic.alt, big.sig = TRUE){
    
    numNulltests = length(statistic.null)
    if(big.sig){
        numSig = sapply(statistic.alt, function(x, statistic.null){ return(sum(statistic.null > x))}, statistic.null = statistic.null)
      }else{
        numSig = sapply(statistic.alt, function(x, statistic.null){ return(sum(statistic.null < x))}, statistic.null = statistic.null)
    }
    numEqual = sapply(statistic.alt, function(x, statistic.null){ return(sum(statistic.null == x))}, statistic.null = statistic.null)
    Uval = runif(length(statistic.alt))
    
    pval.list = (numSig + Uval*(numEqual + 1))/(numNulltests + 1)
    return(pval.list)
}



### 100 window

wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/'
siteSize=2048
treatment='Copper'
null=FALSE
strand='both'
#strand='plus'
#strand='minus'
window.size=100
numSam = 6


##library("DESeq")
##library("DESeq2")

## set up working directory 
setwd(wd.path)



## 60, 30, 20, 10 
filter.cut = 0

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".alt.run/output/DESeq2.600.min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".null.run/output/DESeq2.600.min.pval.", filter.cut, ".txt"))[,1])
 
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))

pval.deseq = get.pval.from.empirical.null.dist.discrete(statistic.null = deseq.null[-del.ix.deseq], statistic.alt = deseq.alt[-del.ix.deseq], big.sig = FALSE)


num.tests = length(deseq.alt)
pval.deseq.600.0 = rep(NA, num.tests)
ix.final = (1:num.tests)[-del.ix.deseq]
pval.deseq.600.0[ix.final] = pval.deseq
length(del.ix.deseq)
## 114617



out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".alt.run/output/DESeq2.pval.discrete.", 600, ".Robj")

save("pval.deseq.600.0", file = out.path)







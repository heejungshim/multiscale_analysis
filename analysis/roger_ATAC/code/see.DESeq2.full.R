#!/usr/bin/env Rscript

## Aim : run DESeq2 with 2048bp windows. I'll focus on `Copper.2048.both'. 
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+


#####################################
## Copy scripts from `run.DESeq.on.roger.ATACseq.R' which contains scrits to run DESeq on Roger's ATAC data.
##################################### 

wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/'
siteSize=2048
treatment='Copper'
null=FALSE
strand='both'
#strand='plus'
#strand='minus'
window.size=100
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


deseq.data.combine = matrix(data=NA, nc = numSam, nr = sum(numSites.list))
for(ss in 1:numSam){
  tmp.dat = matrix(deseq.data[,ss], ncol= numC, byrow =T)
  deseq.data.combine[,ss] = apply(tmp.dat,1,sum)
}




## prepare input format for deseq
condition=factor(c(rep("treated", numSam/2),rep("untreated", numSam/2)))
## prepare data 
deseq.full.data =newCountDataSet(deseq.data.combine,condition)


## let's try 60, 30, 20, 10, 0

#cut.val = c(60, 30, 20, 10, 0)
cut.val = 0

for(cc in 1:length(cut.val)){
  
countData = counts ( deseq.full.data )

###Converts the colData table to a dataframe with vector labels
colData <- data.frame(row.names= c("T1","T2", "T3", "C1", "C2", "C3") ,t=c("T","T","T","C","C","C"),r=as.factor(c(1,2,3,1,2,3)))
### should be a factor!


rsum = rowSums ( countData)
filter.cut = cut.val[cc]
use = (rsum > filter.cut)
countData.filtered = countData[ use, ]

###Perform DESeq2 Analysis

ddsTvC <- DESeqDataSetFromMatrix(
	countData=countData.filtered,
	colData=colData,
	design = ~ r + t)
dds <- DESeq(ddsTvC)
res <- results(dds)


pval.vec = rep(NA, length(use))
pval.vec[use==TRUE] = res$pvalue
pval.filtered = pval.vec
    
## try to save output
## make an output directory
output.dir.path = paste0(input.dir.path, "output/")
if(!file.exists(output.dir.path)){
    dir.create(output.dir.path)
}


## output pval
write.table(pval.filtered, file = paste0(output.dir.path, "DESeq2.full.pval.", filter.cut, ".txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)


}








################################
################################
## Repeat for null
################################
################################


wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/'
siteSize=2048
treatment='Copper'
null=TRUE
strand='both'
#strand='plus'
#strand='minus'
window.size=100
numSam = 6


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


deseq.data.combine = matrix(data=NA, nc = numSam, nr = sum(numSites.list))
for(ss in 1:numSam){
  tmp.dat = matrix(deseq.data[,ss], ncol= numC, byrow =T)
  deseq.data.combine[,ss] = apply(tmp.dat,1,sum)
}




## prepare input format for deseq
condition=factor(c(rep("treated", numSam/2),rep("untreated", numSam/2)))
## prepare data 
deseq.full.data =newCountDataSet(deseq.data.combine,condition)


## let's try 60, 30, 20, 10, 0

#cut.val = c(60, 30, 20, 10, 0)
cut.val = 0

for(cc in 1:length(cut.val)){
  
countData = counts ( deseq.full.data )

###Converts the colData table to a dataframe with vector labels
colData <- data.frame(row.names= c("T1","T2", "T3", "C1", "C2", "C3") ,t=c("T","T","T","C","C","C"),r=as.factor(c(1,2,3,1,2,3)))
### should be a factor!


rsum = rowSums ( countData)
filter.cut = cut.val[cc]
use = (rsum > filter.cut)
countData.filtered = countData[ use, ]

###Perform DESeq2 Analysis

ddsTvC <- DESeqDataSetFromMatrix(
	countData=countData.filtered,
	colData=colData,
	design = ~ r + t)
dds <- DESeq(ddsTvC)
res <- results(dds)


pval.vec = rep(NA, length(use))
pval.vec[use==TRUE] = res$pvalue
pval.filtered = pval.vec
    
## try to save output
## make an output directory
output.dir.path = paste0(input.dir.path, "output/")
if(!file.exists(output.dir.path)){
    dir.create(output.dir.path)
}


## output pval
write.table(pval.filtered, file = paste0(output.dir.path, "DESeq2.full.pval.", filter.cut, ".txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)


}




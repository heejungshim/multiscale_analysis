## `run.DESeq2.logFC.N.LR.on.roger.ATACseq.R' contains scrits to run DESeq2 on Roger's ATAC data and output logFC and LR.
## 
##
## Example Usage (see command in /mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/Copper.2048.plus.100.alt.run/com/deseq.sh) : R CMD BATCH --no-save --no-restore "--args wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/' siteSize=2048 treatment='Copper' null=FALSE strand='plus' window.size=1024 numSam=6 filter.cut=0 test.option=c('Wald', 'LRwithoutS', 'LRwithS')" /mnt/lustre/home/shim/multiscale_analysis/src/R/run.DESeq2.logFC.N.LR.on.roger.ATACseq.R
##
##
## wd.path : working directory path
## siteSize : site size
## treatment : treatment name
## null : indicate whether it's null (control 1 vs control 2) or alternative data
## strand : 'both', 'plus', 'minus'; add two strands, use + strand, or use - strand
## window.size : window size we consider for DESeq analysis
## numSam : number of samples
## filter.cut : analysis includes window with read count > filter.cut
## test.option : what test option to use? c('Wald', 'LRwithoutS', 'LRwithS')
##
##
## Copyright (C) 2014 Heejung Shim
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.


##wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/deseq/'
##siteSize=1024
##treatment='Copper'
##null=TRUE
##strand='both'
##strand='plus'
##strand='minus'
##window.size=1024
##numSam = 6
##filter.cut = 0
##test.option=c('Wald', 'LRwithS', 'LRwithoutS')

args = (commandArgs(TRUE))
eval(parse(text=args))

## check if test.option is correct. 
possible.test.option=c('Wald', 'LRwithoutS', 'LRwithS')
in.test.option = intersect(test.option, possible.test.option)
if(length(in.test.option) == 0){
  cat("Element in test option should be one of these: 'Wald', 'LRwithoutS', 'LRwithS'")
}else{

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
path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/locus/", treatment, ".", siteSize, ".numSites.txt")
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

 
## Converts the colData table to a dataframe with vector labels
colData <- data.frame(row.names= c("T1","T2", "T3", "C1", "C2", "C3") ,t=c("T","T","T","C","C","C"),r=as.factor(c(1,2,3,1,2,3)))

## filter data
rsum = rowSums ( deseq.data)
use = (rsum > filter.cut)
countData.filtered = deseq.data[ use, ]


if(length(intersect('Wald', in.test.option)) > 0){

## Perform DESeq2 Analysis with Wald test option 
ddsTvC <- DESeqDataSetFromMatrix(
	countData=countData.filtered,
	colData=colData,
	design = ~ r + t)
dds <- DESeq(ddsTvC)
res <- results(dds)

mlogLR = qchisq(1-res@listData$pvalue, 1)/2
logFC = log(2)*res@listData$log2FoldChange
selogFC = log(2)*res@listData$lfcSE

mlogLR.vec = logFC.vec = selogFC.vec = rep(NA, length(use))
mlogLR.vec[use==TRUE] = mlogLR
logFC.vec[use==TRUE] = logFC
selogFC.vec[use==TRUE] = selogFC

## try to save output
## make an output directory
if(filter.cut == 0){
  output.dir.path = paste0(input.dir.path, "output/")
  if(!file.exists(output.dir.path)){
    dir.create(output.dir.path)
  }
}else{
  output.dir.path = paste0(input.dir.path, "output", filter.cut, "/")
  if(!file.exists(output.dir.path)){
    dir.create(output.dir.path)
  }
}
  
write.table(mlogLR.vec, file = paste0(output.dir.path, "/logLR.Wald.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)
write.table(logFC.vec, file = paste0(output.dir.path, "/logFC.Wald.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)
write.table(selogFC.vec, file = paste0(output.dir.path, "/selogFC.Wald.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)

## output object just in case we want to take a close look!
save("res", "use", file =paste0(output.dir.path, "/res.Wald.Robj"))

}

if(length(intersect('LRwithoutS', in.test.option)) > 0){

## Perform DESeq2 Analysis with LR without shrinkage
ddsTvC <- DESeqDataSetFromMatrix(
	countData=countData.filtered,
	colData=colData,
	design = ~ r + t)

dds <- DESeq(ddsTvC, test = "LRT", full = design(ddsTvC), reduced = ~ r) 
res <- results(dds)

mlogLR = res@listData$stat/2
logFC = log(2)*res@listData$log2FoldChange
selogFC = log(2)*res@listData$lfcSE

mlogLR.vec = logFC.vec = selogFC.vec = rep(NA, length(use))
mlogLR.vec[use==TRUE] = mlogLR
logFC.vec[use==TRUE] = logFC
selogFC.vec[use==TRUE] = selogFC

## try to save output
## make an output directory
if(filter.cut == 0){
  output.dir.path = paste0(input.dir.path, "output/")
  if(!file.exists(output.dir.path)){
    dir.create(output.dir.path)
  }
}else{
  output.dir.path = paste0(input.dir.path, "output", filter.cut, "/")
  if(!file.exists(output.dir.path)){
    dir.create(output.dir.path)
  }
}
  
write.table(mlogLR.vec, file = paste0(output.dir.path, "/logLR.LRwithoutS.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)
write.table(logFC.vec, file = paste0(output.dir.path, "/logFC.LRwithoutS.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)
write.table(selogFC.vec, file = paste0(output.dir.path, "/selogFC.LRwithoutS.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)

## output object just in case we want to take a close look!
save("res", "use", file =paste0(output.dir.path, "/res.LRwithoutS.Robj"))

}

if(length(intersect('LRwithS', in.test.option)) > 0){

## Perform DESeq2 Analysis with LR with shrinkage
ddsTvC <- DESeqDataSetFromMatrix(
	countData=countData.filtered,
	colData=colData,
	design = ~ r + t)

dds <- DESeq(ddsTvC, test = "LRT", betaPrior = TRUE, full = design(ddsTvC), reduced = ~ r) 
res <- results(dds)

mlogLR = res@listData$stat/2
logFC = log(2)*res@listData$log2FoldChange
selogFC = log(2)*res@listData$lfcSE

mlogLR.vec = logFC.vec = selogFC.vec = rep(NA, length(use))
mlogLR.vec[use==TRUE] = mlogLR
logFC.vec[use==TRUE] = logFC
selogFC.vec[use==TRUE] = selogFC

## try to save output
## make an output directory
if(filter.cut == 0){
  output.dir.path = paste0(input.dir.path, "output/")
  if(!file.exists(output.dir.path)){
    dir.create(output.dir.path)
  }
}else{
  output.dir.path = paste0(input.dir.path, "output", filter.cut, "/")
  if(!file.exists(output.dir.path)){
    dir.create(output.dir.path)
  }
}
  
write.table(mlogLR.vec, file = paste0(output.dir.path, "/logLR.LRwithS.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)
write.table(logFC.vec, file = paste0(output.dir.path, "/logFC.LRwithS.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)
write.table(selogFC.vec, file = paste0(output.dir.path, "/selogFC.LRwithS.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)

## output object just in case we want to take a close look!
save("res", "use", file =paste0(output.dir.path, "/res.LRwithS.Robj"))

}

}

## `run.DESeq2.on.roger.ATACseq.R' contains scrits to run DESeq2 on Roger's ATAC data.
## 
##
## Example Usage (see command in /mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/Copper.2048.plus.100.alt.run/com/deseq.sh) : /data/tools/R-3.1.0/bin/R CMD BATCH --no-save --no-restore "--args wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/' siteSize=2048 treatment='Copper' null=FALSE strand='plus' window.size=100 numSam = 6 filter.cut=0" /mnt/lustre/home/shim/multiscale_analysis/src/R/run.DESeq2.on.roger.ATACseq.R
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
##null=FALSE
##strand='both'
##strand='plus'
##strand='minus'
##window.size=100
##numSam = 6
##filter.cut = 0

args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))
eval(parse(text=args[[2]]))
eval(parse(text=args[[3]]))
eval(parse(text=args[[4]]))
eval(parse(text=args[[5]]))
eval(parse(text=args[[6]]))
eval(parse(text=args[[7]]))
eval(parse(text=args[[8]]))

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

## Perform DESeq2 Analysis
ddsTvC <- DESeqDataSetFromMatrix(
	countData=countData.filtered,
	colData=colData,
	design = ~ r + t)
dds <- DESeq(ddsTvC)
res <- results(dds)

## get p-value
pval.vec = rep(NA, length(use))
pval.vec[use==TRUE] = res$pvalue
pval.filtered =matrix(pval.vec,ncol=numC,byrow=T)

## get minimum p-value for each site
min.pval=apply(pval.filtered,1,min,na.rm=TRUE)
min.pval[is.infinite(min.pval)] = NA

## try to save output
## make an output directory
output.dir.path = paste0(input.dir.path, "output/")
if(!file.exists(output.dir.path)){
    dir.create(output.dir.path)
}

## output min.pval
write.table(min.pval, file = paste0(output.dir.path, "/min.pval.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)

## output results from error check
write.table(sum(err.check==FALSE), file = paste0(output.dir.path, "/errcheck.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)
write.table(sum(is.na(deseq.data)), file = paste0(output.dir.path, "/errcheck.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
#sum(err.check==FALSE)
# 22
#sum(is.na(deseq.data))
# 0

## output object just in case we want to take a close look!
save("res", "use", file =paste0(output.dir.path, "/res.Robj"))















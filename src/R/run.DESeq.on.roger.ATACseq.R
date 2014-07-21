## `run.DESeq.on.roger.ATACseq.R' contains scrits to run DESeq on Roger's ATAC data.
## 
##
## Example Usage (see command in /mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/Copper.2048.plus.100.alt.run/com/deseq.sh) : /data/tools/R-3.1.0/bin/R CMD BATCH --no-save --no-restore "--args wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/' siteSize=2048 treatment='Copper' null=FALSE strand='plus' window.size=100 numSam = 6" /mnt/lustre/home/shim/multiscale_analysis/src/R/run.DESeq.on.roger.ATACseq.R
##
##
## wd.path : working directory path
## siteSize : site size
## treatment : treatment name
## null : indicate whether it's null (control 1 vs control 2) or alternative data
## strand : 'both', 'plus', 'minus'; add two strands, use + strand, or use - strand
## window.size : window size we consider for DESeq analysis
## numSam : number of samples
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



#wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/'
#siteSize=2048
#treatment='Copper'
#null=FALSE
#strand='both'
#strand='plus'
#strand='minus'
#window.size=100
#numSam = 6


args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))
eval(parse(text=args[[2]]))
eval(parse(text=args[[3]]))
eval(parse(text=args[[4]]))
eval(parse(text=args[[5]]))
eval(parse(text=args[[6]]))
eval(parse(text=args[[7]]))

library("DESeq")

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

## prepare input format for deseq
condition=factor(c(rep("treated", numSam/2),rep("untreated", numSam/2)))
## prepare data 
deseq.full.data =newCountDataSet(deseq.data,condition)

## estimate size factors
deseq.full.data=estimateSizeFactors(deseq.full.data)
#sizeFactors(deseq.full.data)
#        1         2         3         4         5         6 
#0.8735805 0.8908987 0.9346553 1.0000000 0.9634925 1.0000000 

## estimate dispersion parameters
deseq.full.data=estimateDispersions(deseq.full.data)

## perform test
resDESeq=nbinomTest(deseq.full.data,"treated", "untreated")
pval=matrix(resDESeq$pval,ncol=numC,byrow=T)

## get minimum p-value for each site
min.pval=apply(pval,1,min,na.rm=TRUE)

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
save("resDESeq", file =paste0(output.dir.path, "/res.Robj"))















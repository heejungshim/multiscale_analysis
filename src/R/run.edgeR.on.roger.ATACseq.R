## `run.edgeR.on.roger.ATACseq.R' contains scrits to run edgeR on Roger's ATAC data.
## 
##
## Example Usage (see command in /mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/edgeR/Copper.2048.plus.100.alt.run/com/edgeR.sh) : /data/tools/R-3.1.0/bin/R CMD BATCH --no-save --no-restore "--args wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/edgeR/' siteSize=2048 treatment='Copper' null=FALSE strand='plus' window.size=100 numSam = 6 filter.cut=0" /mnt/lustre/home/shim/multiscale_analysis/src/R/run.edgeR.on.roger.ATACseq.R
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


##wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/edgeR/'
##siteSize=1024
##treatment='Copper'
##null=FALSE
##strand='both'
##strand='plus'
##strand='minus'
##window.size=300
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

library("edgeR")

## Getting library read depth 
## Path to library read depth
library.read.depth.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/info/"
## getting sample names

## assigen treatment and control name according to input
## treatment    alt     null    control 
## Copper       N702    N705    N706
## Selenium     N703    N705    N706
## Retinoic     N704    N706    N705
name.treatment = NULL
name.control = NULL
if(treatment=='Copper'){
    name.control = "N706"
    if(!null){
        name.treatment = "N702"
    }else{
        name.treatment = "N705"
    }
}
if(treatment=='Selenium'){
    name.control = "N706"
    if(!null){
        name.treatment = "N703"
    }else{
        name.treatment = "N705"
    }
}
if(treatment=='Retinoic'){
    name.control = "N705"
    if(!null){
        name.treatment = "N704"
    }else{
        name.treatment = "N706"
    }
}
## sample names
names.Sam = c("N501", "N502", "N503")
sample.names = c(paste0(name.treatment, names.Sam), paste0(name.control, names.Sam))
## read library read depth 
library.read.depth = rep(0, numSam)
## for reverse
if((strand=='both') | (strand=='plus')){
  path.read.depth = paste0(library.read.depth.path, "library.read.depth.fwd")
  library.read.depth.dat = read.table(path.read.depth, as.is=TRUE)
  for(i in 1:numSam){
    library.read.depth[i] = library.read.depth[i] + library.read.depth.dat[which(library.read.depth.dat[,1] == sample.names[i]),2]
  }
}
## for forward
if((strand=='both') | (strand=='minus')){
  path.read.depth = paste0(library.read.depth.path, "library.read.depth.rev")
  library.read.depth.dat = read.table(path.read.depth, as.is=TRUE)
  for(i in 1:numSam){
    library.read.depth[i] = library.read.depth[i] + library.read.depth.dat[which(library.read.depth.dat[,1] == sample.names[i]),2]
  }
}


## set up working directory 
setwd(wd.path)

## input data path (input data are under deseq directory
temp = unlist(strsplit(wd.path, "/"))
input.path = paste(c(temp[-length(temp)], "deseq/"), collapse = "/")

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
input.dir.path = paste0(input.path, dir.name, ".run/")
numC = siteSize%/%window.size
edgeR.data = matrix(data=NA, nc = numSam, nr = sum(numSites.list)*numC)
st.chr.ix = NULL
en.chr.ix = NULL
st.chr.ix[1] = 1
err.check = rep(NA, 22)
for(chr in 1:22){
    #chr = 1
    numSites = numSites.list[chr]
    numRow = numSites*numC
    en.chr.ix[chr] = st.chr.ix[chr] + numRow - 1
    edgeR.data.each = read.table(paste0(input.dir.path, "data.", chr, ".txt"))
    if(dim(edgeR.data.each)[1] == numRow){
        err.check[chr] = FALSE
        edgeR.data[st.chr.ix[chr]:en.chr.ix[chr],] = as.matrix(edgeR.data.each)
    }
    st.chr.ix[chr+1] = en.chr.ix[chr] + 1
}
#sum(err.check==FALSE)
# 22
#sum(is.na(edgeR.data))
# 0


## filter data
rsum = rowSums(edgeR.data)
use = (rsum > filter.cut)
countData.filtered = edgeR.data[ use, ]

## Perform edgeR Analysis
gp = factor(c(rep(0,3), rep(1,3)) ,labels=c("A","B"))
res = DGEList(counts=countData.filtered, group=gp, lib.size=library.read.depth)
res = calcNormFactors(res, method="RLE")
res = estimateCommonDisp(res)
res = estimateTagwiseDisp(res)
res = exactTest(res, dispersion="auto")

##str(res$table)
##'data.frame':	726299 obs. of  3 variables:
## $ logFC : num  0.3893 0.5911 -1.8901 -0.0487 0.4558 ...
## $ logCPM: num  -2.23 -1.7 -5.06 -4.29 -3.08 ...
## $ PValue: num  0.35 0.119 1 1 0.467 ...

## get p-value
pval.vec = rep(NA, length(use))
pval.vec[use==TRUE] = res$table$PValue
pval.filtered =matrix(pval.vec,ncol=numC,byrow=T)

## get minimum p-value for each site
min.pval=apply(pval.filtered,1,min,na.rm=TRUE)
min.pval[is.infinite(min.pval)] = NA

## try to save output
## make an output directory
output.dir.path = paste0(wd.path, dir.name, ".run/")
if(!file.exists(output.dir.path)){
    dir.create(output.dir.path)
}
output.dir.path = paste0(output.dir.path, "output/")
if(!file.exists(output.dir.path)){
    dir.create(output.dir.path)
}


## output min.pval
write.table(min.pval, file = paste0(output.dir.path, "/min.pval.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)

## output results from error check
write.table(sum(err.check==FALSE), file = paste0(output.dir.path, "/errcheck.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)
write.table(sum(is.na(edgeR.data)), file = paste0(output.dir.path, "/errcheck.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
#sum(err.check==FALSE)
# 22
#sum(is.na(deseq.data))
# 0

## output object just in case we want to take a close look!
save("res", "use", file =paste0(output.dir.path, "/res.Robj"))

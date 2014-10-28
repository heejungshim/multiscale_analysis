#!/usr/bin/env Rscript

## Aim : Resutls from DESeq analysis are very bad. I'm not sure I used DESseq properly. So I'll try different ways to use DESeq and see if results make sens. For this purpose, I'll focus on `Copper.2048.both.300'. I'll try to exclude 300bp windows with low read counts.
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
window.size=300
numSam = 6


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

##############################
## select filtering criteria
## 0, 5, 10, 20, 30, 60
##############################
rsum = rowSums ( counts ( deseq.full.data ))
filter.cut = 60
use = (rsum > filter.cut)
all.use = apply(matrix(use,ncol=numC,byrow=T),1,sum)
sum(all.use > 0)
length(all.use)
sum(all.use > 0)/length(all.use)

## filter.cut : 0
## 186522
## 197969
## 0.9421778

## filter.cut : 5
## 185843
## 197969
## 0.938748

## filter.cut : 10
## 185328
## 197969
## 0.9361466

## filter.cut : 20
## 184034
## 197969
## 0.9296102

## filter.cut : 30
## 179342
## 197969
## 0.9059095

## filter.cut : 60
## 111065
## 197969
## 0.5610222

## let's try 60, 30, 20, 10

rsum = rowSums ( counts ( deseq.full.data ))
filter.cut = 10
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
write.table(min.pval, file = paste0(output.dir.path, "/min.pval.", filter.cut, ".txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)


#############################
## Without filtering
#############################

deseq.full.data = estimateSizeFactors(deseq.full.data)
sizeFactors(deseq.full.data)

## estimate dispersion parameters
deseq.full.data=estimateDispersions(deseq.full.data)

## perform test
resDESeq=nbinomTest(deseq.full.data,"treated", "untreated")
pval=matrix(resDESeq$pval,ncol=numC,byrow=T)

## get minimum p-value for each site
min.pval=apply(pval,1,min,na.rm=TRUE)
min.pval[is.infinite(min.pval)] = NA

## try to save output
## make an output directory
output.dir.path = paste0(input.dir.path, "output/")
if(!file.exists(output.dir.path)){
    dir.create(output.dir.path)
}

## output min.pval
write.table(min.pval, file = paste0(output.dir.path, "/min.pval.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)



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
window.size=300
numSam = 6


## library("DESeq")

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

##############################
## select filtering criteria
## 0, 5, 10, 20, 30, 60
##############################
rsum = rowSums ( counts ( deseq.full.data ))
filter.cut = 60
use = (rsum > filter.cut)
all.use = apply(matrix(use,ncol=numC,byrow=T),1,sum)
sum(all.use > 0)
length(all.use)
sum(all.use > 0)/length(all.use)

## filter.cut : 0
## 186507
## 197969
## 0.942102

## filter.cut : 5
## 185738
## 197969
## 0.9382176

## filter.cut : 10
## 185137
## 197969
## 0.9351818

## filter.cut : 20
## 182506
## 197969
## 0.9218918

## filter.cut : 30
## 172008
## 197969
## 0.8688633

## filter.cut : 60
## 90598
## 197969
## 0.4576373


## let's try 60, 30, 20, 10

rsum = rowSums ( counts ( deseq.full.data ))
filter.cut = 10
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
write.table(min.pval, file = paste0(output.dir.path, "/min.pval.", filter.cut, ".txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)


#############################
## Without filtering
#############################

deseq.full.data = estimateSizeFactors(deseq.full.data)
sizeFactors(deseq.full.data)

## estimate dispersion parameters
deseq.full.data=estimateDispersions(deseq.full.data)

## perform test
resDESeq=nbinomTest(deseq.full.data,"treated", "untreated")
pval=matrix(resDESeq$pval,ncol=numC,byrow=T)

## get minimum p-value for each site
min.pval=apply(pval,1,min,na.rm=TRUE)
min.pval[is.infinite(min.pval)] = NA

## try to save output
## make an output directory
output.dir.path = paste0(input.dir.path, "output/")
if(!file.exists(output.dir.path)){
    dir.create(output.dir.path)
}

## output min.pval
write.table(min.pval, file = paste0(output.dir.path, "/min.pval.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)





###################################################
## let's compare results from three methods
## copy scripts from see.statistic.3method.both.R
#####################################################

#######################################
## Step 1: Get p-value from ms and wave
#######################################

all.name = paste0(treatment,".", siteSize, ".", strand)

deseq.300.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.txt"))[,1])
length(deseq.300.alt)

deseq.300.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.txt"))[,1])
length(deseq.300.null)

#> which(is.na(ms.alt)==TRUE)
#[1] 35550
#> which(is.na(ms.null)==TRUE)
#[1]  61034 126327 126328 133432
del.ix = union(union(which(is.na(deseq.300.alt)==TRUE), which(is.na(deseq.300.null)==TRUE)),  c(35550, 61034, 126327, 126328, 133432))


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/wave.new.", all.name, ".Robj"))
load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/ms.new.", all.name, ".Robj"))
load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/deseq.", all.name, ".Robj"))



## length(pval.ms.new)
## 186399
## length(pval.wave.new)
## 186399
## length(del.ix)
## 11570
## length(deseq.300.null)
## 197969
## 197969 - 11570
## 186399
## length(pval.deseq.300)
## 186399

num.tests = length(deseq.300.null)
pval.ms = rep(NA, num.tests)
pval.wave = rep(NA, num.tests)
pval.deseq.3 = rep(NA, num.tests)
ix.final = (1:num.tests)[-del.ix]
length(ix.final)
pval.ms[ix.final] = pval.ms.new
pval.wave[ix.final] = pval.wave.new
pval.deseq.3[ix.final] = pval.deseq.300



out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/pval.ms.wave.", all.name, ".Robj")
save("pval.ms", "pval.wave", file = out.path)




##################################################################
# Step 2: get p-value from empirical null distribution of test statistic (DESeq)
##################################################################

get.pval.from.empirical.null.dist <- function(statistic.null, statistic.alt, big.sig = TRUE){
    
    numNulltests = length(statistic.null)
    if(big.sig){
        numSig = sapply(statistic.alt, function(x, statistic.null){ return(sum(statistic.null >= x))}, statistic.null = statistic.null)
    }else{
        numSig = sapply(statistic.alt, function(x, statistic.null){ return(sum(statistic.null <= x))}, statistic.null = statistic.null)
    }
    pval.list = sapply(numSig, function(x, numNulltests){ return(runif(1, (x+1)/(numNulltests+2), (x+1)/(numNulltests+1)))}, numNulltests = numNulltests)

    return(pval.list)
}


## 60, 30, 20, 10 
filter.cut = 60

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.", filter.cut, ".txt"))[,1])
 
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))

pval.deseq = get.pval.from.empirical.null.dist(statistic.null = deseq.null[-del.ix.deseq], statistic.alt = deseq.alt[-del.ix.deseq], big.sig = FALSE)

pval.deseq.3.60 = rep(NA, num.tests)
ix.final = (1:num.tests)[-del.ix.deseq]
pval.deseq.3.60[ix.final] = pval.deseq
length(del.ix.deseq)
## 114617




filter.cut = 30

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.", filter.cut, ".txt"))[,1])
 
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))

pval.deseq = get.pval.from.empirical.null.dist(statistic.null = deseq.null[-del.ix.deseq], statistic.alt = deseq.alt[-del.ix.deseq], big.sig = FALSE)

pval.deseq.3.30 = rep(NA, num.tests)
ix.final = (1:num.tests)[-del.ix.deseq]
pval.deseq.3.30[ix.final] = pval.deseq
length(del.ix.deseq)
## 27443



filter.cut = 20

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.", filter.cut, ".txt"))[,1])
 
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))

pval.deseq = get.pval.from.empirical.null.dist(statistic.null = deseq.null[-del.ix.deseq], statistic.alt = deseq.alt[-del.ix.deseq], big.sig = FALSE)

pval.deseq.3.20 = rep(NA, num.tests)
ix.final = (1:num.tests)[-del.ix.deseq]
pval.deseq.3.20[ix.final] = pval.deseq
length(del.ix.deseq)
## 15831



filter.cut = 10

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.", filter.cut, ".txt"))[,1])
 
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))

pval.deseq = get.pval.from.empirical.null.dist(statistic.null = deseq.null[-del.ix.deseq], statistic.alt = deseq.alt[-del.ix.deseq], big.sig = FALSE)

pval.deseq.3.10 = rep(NA, num.tests)
ix.final = (1:num.tests)[-del.ix.deseq]
pval.deseq.3.10[ix.final] = pval.deseq
length(del.ix.deseq)
## 13014


out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/deseq.pval.Robj")

save("pval.deseq.3", "pval.deseq.3.60", "pval.deseq.3.30", "pval.deseq.3.20", "pval.deseq.3.10", file = out.path)












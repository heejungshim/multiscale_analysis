## `get.effect.on.roger.newATACseq.R' compute effect size from multiseq for selected sites
#
## Example Usage : R CMD BATCH --no-save --no-restore "--args info.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/summary/Copper.1024.both.msOnly.ms.DESeq300.info' DESeq.100.info.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/deseq/Copper.1024.both.100.alt.run/output/res.Robj' DESeq.300.info.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/deseq/Copper.1024.both.300.alt.run/output/res.Robj' DESeq.full.info.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/deseq/Copper.1024.both.1024.alt.run/output/res.Robj' out.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/code/' wave.out.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/wave/' file.name='ES' siteSize=1024 treatment='Copper' strand='both' sig.level=2 wave.effect=TRUE multiseq.effect=TRUE deseq.100.effect=TRUE deseq.300.effect=TRUE" /mnt/lustre/home/shim/multiscale_analysis/src/R/make.effect.figure.on.roger.ATACseq.R
##
##
## info.path : path to file that contains information on sites of interest (index chr sites st.posi en.posi, ...)
## pcr.path : path to directory where position with prc artifacts are saved
## out.path : path to directory where effect size will be saved
## siteSize : site size
## treatment : treatment name
## strand : 'both', 'plus', 'minus'; add two strands, use + strand, or use - strand
##
## Copyright (C) 2015 Heejung Shim
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


setwd("/mnt/lustre/home/shim/multiscale_analysis")

library("multiseq")
library("ashr")

multiscale.analysis.repodir <- scan(".multiscale_analysis.repodir.txt", what=character())
source(paste0(multiscale.analysis.repodir, "/src/R/utils.R"))
source(paste0(multiscale.analysis.repodir, "/src/R/my.utils.R"))

info.path = '/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/summary/EffectSize/Copper.1024.both.f11.sites'
pcr.path = '/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/multiscale/Copper.1024.both.alt.output/'
out.path = '/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/summary/EffectSize/Copper.1024.both/'
siteSize=1024
treatment='Copper'
strand='both'


args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))
eval(parse(text=args[[2]]))
eval(parse(text=args[[3]]))
eval(parse(text=args[[4]]))
eval(parse(text=args[[5]]))
eval(parse(text=args[[6]]))


## assigen treatment and control name according to input
## treatment    alt     null    control 
## Copper       N702    N705    N706
## Selenium     N703    N705    N706
## Retinoic     N704    N706    N705

##############################################
## sample name and sample file for alternative data
##############################################
null = FALSE
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

## Make a list of sample names and a list of hdf5 file names : treatment first and control later.
sample.names = c(paste0(name.treatment, names.Sam), paste0(name.control, names.Sam))
sample.files = paste0(sample.names, ".qfiltered10")

sample.names.alt = sample.names
sample.files.alt = sample.files


## Path to directory which contain ATAC-seq data as hdf5 format, 
hdf5.data.path = "/data/share/genome_db/hg19/roger_atacseq2/"

## Make a covariate
g = c(rep(0, length(names.Sam)), rep(1, length(names.Sam)))

## Path to library read depth
library.read.depth.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/info/"

## set up working directory
setwd(out.path)

## read information on seletec regions
dat.info = read.table(file=info.path, header = TRUE, as.is=TRUE)

numSites = dim(dat.info)[1]

for(ss in 1:numSites){

  ss = 8
## read location information
chr = dat.info$chr[ss]
st.posi = dat.info$st.posi[ss]
en.posi = dat.info$en.posi[ss]
site = dat.info$sites[ss]

## read pcr posi
pcr.posi = vector("list", 2)
pcr.file.path = paste0(pcr.path, "pcrposi.", chr, ".", site, ".out")
if(file.exists(pcr.file.path)){
  pcr.posi[[1]] = as.numeric(read.table(file = pcr.file.path, nrows = 1))[-1]
  pcr.posi[[2]] = as.numeric(read.table(file = pcr.file.path, nrows = 1, skip = 1))[-1]
}

#####################
## read data alt
#####################
null = FALSE
sample.names = sample.names.alt
sample.files = sample.files.alt
numSam = length(sample.names)
numBPs = siteSize
library.read.depth = rep(0, numSam)
ATAC.dat = matrix(data=0, nr = numSam, nc = numBPs)
ATAC.dat.orig = matrix(data=0, nr = numSam, nc = numBPs)
ATAC.dat.full = matrix(data=0, nr = numSam, nc = numBPs)

pcr.ix = 1

## for fwd
if((strand=='both') | (strand=='plus')){

    ## read library read depth
    path.read.depth = paste0(library.read.depth.path, "library.read.depth.fwd")
    library.read.depth.dat = read.table(path.read.depth, as.is=TRUE)
    for(i in 1:numSam){
        library.read.depth[i] = library.read.depth[i] + library.read.depth.dat[which(library.read.depth.dat[,1] == sample.names[i]),2]
    }

    ## read read counts for a given region
    ## for + strand, we need get reads at locations that are shifted 4bp to left
    ATAC.dat.fwd = matrix(data=NA, nr = numSam, nc = numBPs)
    for(i in 1:numSam){
        path.fwd = paste0(hdf5.data.path, sample.files[i] , ".fwd.h5")
        ATAC.dat.fwd[i, 1:numBPs] = as.matrix(get.counts.h5(path.fwd, paste0("chr", chr), st.posi-4, en.posi-4))
    }

    ## remove pcr artifacts
    pcr.removed.fwd = remove.pcr.artifacts.in.known.posi(data=ATAC.dat.fwd, known.pcr.posi = pcr.posi[[pcr.ix]] , win.half.size = 50, prop.thresh = 0.9)
    ATAC.dat = ATAC.dat + pcr.removed.fwd$data
    pcr.ix = pcr.ix + 1

    ## should delete
    ATAC.dat.orig = ATAC.dat.orig + ATAC.dat.fwd
    pcr.removed.fwd = remove.pcr.artifacts(data=ATAC.dat.fwd, win.half.size = 50, prop.thresh = 0.9)
    ATAC.dat.full = ATAC.dat.full + pcr.removed.fwd$data
}

## for reverse
if((strand=='both') | (strand=='minus')){

    ## read library read depth
    path.read.depth = paste0(library.read.depth.path, "library.read.depth.rev")
    library.read.depth.dat = read.table(path.read.depth, as.is=TRUE)
    for(i in 1:6){
        library.read.depth[i] = library.read.depth[i] + library.read.depth.dat[which(library.read.depth.dat[,1] == sample.names[i]),2]
    }
    
    ## read read counts for a given region
    ## for - strand, we need get reads at locations that are shifted 4bp to right
    ATAC.dat.rev = matrix(data=NA, nr = numSam, nc = numBPs)
    for(i in 1:numSam){
        path.rev = paste0(hdf5.data.path, sample.files[i] , ".rev.h5")
        ATAC.dat.rev[i, 1:numBPs] = as.matrix(get.counts.h5(path.rev, paste0("chr", chr), st.posi+4, en.posi+4))
    }

    ## remove pcr artifacts
    pcr.removed.rev = remove.pcr.artifacts.in.known.posi(data=ATAC.dat.rev, known.pcr.posi = pcr.posi[[pcr.ix]] , win.half.size=50, prop.thresh=0.9)
    ATAC.dat = ATAC.dat + pcr.removed.rev$data
    pcr.ix = pcr.ix + 1

    ## should delete
    ATAC.dat.orig = ATAC.dat.orig + ATAC.dat.rev
    pcr.removed.rev = remove.pcr.artifacts(data=ATAC.dat.rev, win.half.size = 50, prop.thresh = 0.9)
    ATAC.dat.full = ATAC.dat.full + pcr.removed.rev$data

    
}

sum(ATAC.dat.full != ATAC.dat)
sum(ATAC.dat.full == ATAC.dat)
sum(ATAC.dat.orig == ATAC.dat)
sum(ATAC.dat.orig != ATAC.dat)


  
phenoD = ATAC.dat

## get effect size
genoD = g
res = multiseq(x = phenoD, g = genoD, read.depth = library.read.depth)

effect.mean = -res$effect.mean
effect.sd = sqrt(res$effect.var)
ES.dat = data.frame(effect.mean, effect.sd)

out.each.path = paste0(out.path, "chr", chr, ".", st.posi, ".", en.posi) 
write.table(ES.dat, file = out.each.path, quote=F, row.names = F, col.names = F)

}

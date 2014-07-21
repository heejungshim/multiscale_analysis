## `preprocess.wave.DESeq.N.run.multiscale.on.roger.ATACseq.R' combines `preprocess.wave.DESeq.roger.ATACseq.R' and `run.multiscale.on.roger.ATACseq.R'
##
##
## Example Usage (see command in /mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/gen.data/com) : R CMD BATCH --no-save --no-restore "--args chr=1 sites.ix=1 wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/' siteSize=2048 treatment='Copper' null=FALSE strand='both' meanR.thresh=1 window.size.list=c(100,300) wavelet.preprocess=TRUE deseq.preprocess=TRUE" /mnt/lustre/home/shim/multiscale_analysis/src/R/preprocess.wave.DESeq.N.run.multiscale.on.roger.ATACseq.R
##
##
## chr : chromosome
## sites.ix : default=NULL; if it is null, run multiseq on all sites or we can specifiy partifular site
## wd.path : working directory path
## siteSize : site size
## treatment : treatment name
## null : indicate whether it's null (control 1 vs control 2) or alternative data
## strand : 'both', 'plus', 'minus'; add two strands, use + strand, or use - strand
## meanR.thresh=1 : for wavelet preprocess
## window.size.list=c(100,300) : list of window sizes we consider for DESeq analysis
## wavelet.preprocess=TRUE
## deseq.preprocess=TRUE
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




setwd("/mnt/lustre/home/shim/multiscale_analysis")

library("multiseq")
library("ashr")


multiscale.analysis.repodir <- scan(".multiscale_analysis.repodir.txt", what=character())
source(paste0(multiscale.analysis.repodir, "/src/R/utils.R"))
source(paste0(multiscale.analysis.repodir, "/src/R/my.utils.R"))
WaveQTL.repodir <- scan(".WaveQTL.repodir.txt", what=character())



#chr=1
#sites.ix=1
#wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/'
#siteSize=2048
#treatment='Selenium'

#null=FALSE
#strand='both'
#strand='plus'
#strand='minus'
#meanR.thresh=1
#window.size.list=c(100,300)
#wavelet.preprocess=TRUE
#deseq.preprocess=TRUE




args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))
eval(parse(text=args[[2]]))
eval(parse(text=args[[3]]))
eval(parse(text=args[[4]]))
eval(parse(text=args[[5]]))
eval(parse(text=args[[6]]))
eval(parse(text=args[[7]]))
eval(parse(text=args[[8]]))
eval(parse(text=args[[9]]))
eval(parse(text=args[[10]]))
eval(parse(text=args[[11]]))


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

## Path to directory which contain ATAC-seq data as hdf5 format, 
hdf5.data.path = "/data/share/genome_db/hg19/roger_atacseq/"

## Make a list of sample names and a list of hdf5 file names : treatment first and control later.
sample.names = c(paste0(name.treatment, names.Sam), paste0(name.control, names.Sam))
sample.files = paste0(sample.names, ".qfiltered10")


## Make a covariate
g = c(rep(0, length(names.Sam)), rep(1, length(names.Sam)))

## Path to library read depth
library.read.depth.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/info/"




## set up working directory 
setwd(wd.path)

## make output directory name and output directory
if(!null){
    output.dir.name=paste0(treatment,".", siteSize, ".", strand,  ".alt")
}else{
    output.dir.name=paste0(treatment,".", siteSize, ".", strand, ".null")
}


multiseq.out.dir.path = paste0(wd.path, "multiscale/", output.dir.name, ".output") 
if(!file.exists(multiseq.out.dir.path)){
    dir.create(multiseq.out.dir.path)
}


    
## make warning message directory
warning.dir.path = paste0(wd.path, "multiscale/", output.dir.name, ".warnings") 
if(!file.exists(warning.dir.path)){
    dir.create(warning.dir.path)
}


## make output directory name and output directory for wavelets
if(wavelet.preprocess){

    wave.out.dir.path = paste0(wd.path, "wave/", output.dir.name, ".data/") 
    if(!file.exists(wave.out.dir.path)){
        dir.create(wave.out.dir.path)
    }
}

## make output directory name and output directory for DESeq
if(deseq.preprocess){
    for(ss in 1:length(window.size.list)){
        window.size = window.size.list[ss]
        if(!null){
            output.dir.name=paste0(treatment,".", siteSize, ".", strand,  ".", window.size, ".alt")
        }else{
            output.dir.name=paste0(treatment,".", siteSize, ".", strand, ".", window.size, ".null")
        }
 
        deseq.out.dir.path = paste0(wd.path, "deseq/", output.dir.name, ".data/") 
        if(!file.exists(deseq.out.dir.path)){
            dir.create(deseq.out.dir.path)
        }
    }
}    




#############################
# read location information and numSites for each chromosome
#############################
path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/", treatment, ".", siteSize, ".chr", chr, ".locus")
loc.info = read.table(path, header=TRUE)

path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/", treatment, ".", siteSize, ".numSites.txt")
numSites = scan(path)[chr]


## whether we run a whole chromosome together or not 
if(is.null(sites.ix)){
    st.sites = 1
    en.sites = numSites
    write.table(NULL, file = paste0(multiseq.out.dir.path, "/res.", chr, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE)
    write.table(NULL, file = paste0(multiseq.out.dir.path, "/pcrposi.", chr, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE)
}else{
    st.sites = sites.ix
    en.sites = sites.ix
}

for(sites in st.sites:en.sites){

# sites = 1    
## create file for warning message
warn.path = paste0(warning.dir.path, "/warnings.", chr, ".", sites, ".txt")
warnings.file <- file(warn.path, open="wt")
sink(warnings.file, type="message")



#############################
# read location information 
#############################
st.posi = loc.info$st.posi[sites]
en.posi = loc.info$en.posi[sites]


#####################
# read data
#####################
numSam = length(sample.names)
numBPs = siteSize
library.read.depth = rep(0, numSam)
ATAC.dat = matrix(data=0, nr = numSam, nc = numBPs)

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

    ATAC.dat = ATAC.dat + ATAC.dat.fwd

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

    ATAC.dat = ATAC.dat + ATAC.dat.rev

}


## remove pcr artifacts
pcr.removed = remove.pcr.artifacts(data=ATAC.dat, win.half.size=50, prop.thresh=0.9)
phenoD = pcr.removed$data
pcr.posi = pcr.removed$posi.with.pcr.artifacts



#########################################
# data preprocessing for wavelet analysis
#########################################

if(wavelet.preprocess){

    source(paste0(WaveQTL.repodir, "/R/WaveQTL_preprocess_funcs.R"))
    res = WaveQTL_preprocess(Data = phenoD, library.read.depth = library.read.depth , Covariates = NULL, meanR.thresh = meanR.thresh)

    filteredWCs = res$filtered.WCs
    norm.WCs = res$WCs

    this.path = paste0(wave.out.dir.path, "WC.", chr, ".", sites, ".txt")
    write.table(norm.WCs, file = this.path, quote= FALSE, row.names = FALSE, col.names = FALSE)
    this.path = paste0(wave.out.dir.path, "use.", chr, ".", sites, ".txt")
    cat(filteredWCs, file = this.path)

}

 
#########################################
# data preprocessing for DESeq
#########################################
if(deseq.preprocess){

    for(ss in 1:length(window.size.list)){
        window.size = window.size.list[ss]
        numC = numBPs%/%window.size
        mat = matrix(data=NA, nc = numSam, nr = numC)
        st = 1
        for(c in 1:(numC-1)){
             en = st + window.size  - 1
            mat[c,] = apply(phenoD[,st:en], 1, sum)
            st = en + 1
        }

        en = numBPs
        mat[numC,] = apply(phenoD[,st:en], 1, sum)
  
        ## save deseq data
        if(!null){
            output.dir.name=paste0(treatment,".", siteSize, ".", strand,  ".", window.size, ".alt")
        }else{
            output.dir.name=paste0(treatment,".", siteSize, ".", strand, ".", window.size, ".null")
        }
 
        deseq.out.dir.path = paste0(wd.path, "deseq/", output.dir.name, ".data/")
        this.path = paste0(deseq.out.dir.path, "data.", chr, ".", sites, ".txt")
        write.table(mat, file = this.path, quote= FALSE, row.names = FALSE, col.names = FALSE)
    }

}



## perform test
genoD = g
res = permutation.logLR(pheno.dat = phenoD, geno.dat = genoD, library.read.depth = library.read.depth, numPerm = NULL, use.default.compute.logLR = TRUE, cxx=TRUE)

## write output
if(is.null(sites.ix)){
    write.table(t(c(sites, res$logLR)), file = paste0(multiseq.out.dir.path, "/res.", chr, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    write.table(t(c(sites, pcr.posi)), file = paste0(multiseq.out.dir.path, "/pcrposi.", chr, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}else{
    write.table(res$logLR, file = paste0(multiseq.out.dir.path, "/res.", chr, ".", sites, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE)
    if(is.null(pcr.posi)){
        write.table(pcr.posi, file = paste0(multiseq.out.dir.path, "/pcr.posi.", chr, ".", sites, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE)
    }else{
        write.table(t(pcr.posi), file = paste0(multiseq.out.dir.path, "/pcr.posi.", chr, ".", sites, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE)        
    }
}


}
















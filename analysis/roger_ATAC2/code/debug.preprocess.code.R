## copy from /mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/gen.data/com/Copper.1024.both.null/all.4.sh

chr=4
sites.ix=56
wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/'
siteSize=1024
treatment='Copper'
null=TRUE
strand='both'
meanR.thresh=2
window.size.list=c(100,300,1024)
wavelet.preprocess.QT=TRUE
wavelet.preprocess.NoQT=FALSE
deseq.preprocess=TRUE


## copy from ~/multiscale_analysis/src/R/preprocess.wave.DESeq.N.run.multiscale.on.roger.ATACseq.R

setwd("/mnt/lustre/home/shim/multiscale_analysis")

library("multiseq")
library("ashr")


multiscale.analysis.repodir <- scan(".multiscale_analysis.repodir.txt", what=character())
source(paste0(multiscale.analysis.repodir, "/src/R/utils.R"))
source(paste0(multiscale.analysis.repodir, "/src/R/my.utils.R"))
WaveQTL.repodir <- scan(".WaveQTL.repodir.txt", what=character())


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
hdf5.data.path = "/data/share/genome_db/hg19/roger_atacseq2/"

## Make a list of sample names and a list of hdf5 file names : treatment first and control later.
sample.names = c(paste0(name.treatment, names.Sam), paste0(name.control, names.Sam))
sample.files = paste0(sample.names, ".qfiltered10")


## Make a covariate
g = c(rep(0, length(names.Sam)), rep(1, length(names.Sam)))

## Path to library read depth
library.read.depth.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/info/"




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
##warning.dir.path = paste0(wd.path, "multiscale/", output.dir.name, ".warnings") 
##if(!file.exists(warning.dir.path)){
##    dir.create(warning.dir.path)
##}


## make output directory name and output directory for wavelets
if(wavelet.preprocess.QT){  
    wave.out.dir.path = paste0(wd.path, "wave/", output.dir.name, ".data/") 
    if(!file.exists(wave.out.dir.path)){
        dir.create(wave.out.dir.path)
    }
}

if(wavelet.preprocess.NoQT){  
    waveNoQT.out.dir.path = paste0(wd.path, "waveNoQT/", output.dir.name, ".data/") 
    if(!file.exists(waveNoQT.out.dir.path)){
        dir.create(waveNoQT.out.dir.path)
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
path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/locus/", treatment, ".", siteSize, ".chr", chr, ".locus")
loc.info = read.table(path, header=TRUE)

path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/locus/", treatment, ".", siteSize, ".numSites.txt")
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

##for(sites in st.sites:en.sites){

sites = 56    
## create file for warning message
##warn.path = paste0(warning.dir.path, "/warnings.", chr, ".", sites, ".txt")
##warnings.file <- file(warn.path, open="wt")
##sink(warnings.file, type="message")



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
pcr.posi = vector("list", 2)
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
    pcr.removed.fwd = remove.pcr.artifacts(data=ATAC.dat.fwd, win.half.size=50, prop.thresh=0.9)
    ATAC.dat = ATAC.dat + pcr.removed.fwd$data
    if(!is.null(pcr.removed.fwd$posi.with.pcr.artifacts)){
        pcr.posi[[pcr.ix]] = pcr.removed.fwd$posi.with.pcr.artifacts
    }
    pcr.ix = pcr.ix + 1

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
    pcr.removed.rev = remove.pcr.artifacts(data=ATAC.dat.rev, win.half.size=50, prop.thresh=0.9)
    ATAC.dat = ATAC.dat + pcr.removed.rev$data
    if(!is.null(pcr.removed.rev$posi.with.pcr.artifacts)){
        pcr.posi[[pcr.ix]] = pcr.removed.rev$posi.with.pcr.artifacts
    }
    pcr.ix = pcr.ix + 1
    
}

phenoD = ATAC.dat

#########################################
# data preprocessing for wavelet analysis
#########################################

if(wavelet.preprocess.QT){

    source(paste0(WaveQTL.repodir, "/R/WaveQTL_preprocess_funcs.R"))
    res = WaveQTL_preprocess(Data = phenoD, library.read.depth = library.read.depth , Covariates = NULL, meanR.thresh = meanR.thresh)

    filteredWCs = res$filtered.WCs
    norm.WCs = res$WCs

    this.path = paste0(wave.out.dir.path, "WC.", chr, ".", sites, ".txt")
    write.table(norm.WCs, file = this.path, quote= FALSE, row.names = FALSE, col.names = FALSE)
    this.path = paste0(wave.out.dir.path, "use.", chr, ".", sites, ".txt")
    cat(filteredWCs, file = this.path)
  }

if(wavelet.preprocess.NoQT){

    source(paste0(WaveQTL.repodir, "/R/WaveQTL_preprocess_funcs.R"))
    res = WaveQTL_preprocess(Data = phenoD, library.read.depth = library.read.depth , Covariates = NULL, meanR.thresh = meanR.thresh, no.QT = TRUE)

    filteredWCs = res$filtered.WCs
    norm.WCs = res$WCs

    this.path = paste0(waveNoQT.out.dir.path, "WC.", chr, ".", sites, ".txt")
    write.table(norm.WCs, file = this.path, quote= FALSE, row.names = FALSE, col.names = FALSE)
    this.path = paste0(waveNoQT.out.dir.path, "use.", chr, ".", sites, ".txt")
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
        if(numC > 1){
          for(c in 1:(numC-1)){
            en = st + window.size  - 1
            mat[c,] = apply(phenoD[,st:en], 1, sum)
            st = en + 1
          }
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
###Error in (-npoint):0 : result would be too long a vector


## write output
if(is.null(sites.ix)){
    write.table(t(c(sites, res$logLR)), file = paste0(multiseq.out.dir.path, "/res.", chr, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    write.table(t(c(sites, pcr.posi[[1]])), file = paste0(multiseq.out.dir.path, "/pcrposi.", chr, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    write.table(t(c(sites, pcr.posi[[2]])), file = paste0(multiseq.out.dir.path, "/pcrposi.", chr, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}else{
    write.table(res$logLR, file = paste0(multiseq.out.dir.path, "/res.", chr, ".", sites, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE)
    for(m in 1:(pcr.ix-1)){
        if(m == 1){
            write.table(t(c(m, pcr.posi[[m]])), file = paste0(multiseq.out.dir.path, "/pcrposi.", chr, ".", sites, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE)
        }else{
            write.table(t(c(m, pcr.posi[[m]])), file = paste0(multiseq.out.dir.path, "/pcrposi.", chr, ".", sites, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
        }
    }
}


}






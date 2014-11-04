## `make.effect.figure.on.roger.ATACseq.R' makes effect size figures from multiseq, wavelets, DESeq for selected sites
#
## Example Usage : R CMD BATCH --no-save --no-restore "--args info.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/tmp/Copper.2048.both.msOnly.info' out.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/code/' wave.out.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/wave/' file.name='msOnly' siteSize=2048 treatment='Copper' null=FALSE strand='both' sig.level=2 wave.effect=TRUE multiseq.effect=TRUE deseq.effect=FALSE" /mnt/lustre/home/shim/multiscale_analysis/src/R/make.effect.figure.on.roger.ATACseq.R
##
##
## info.path : path to file that contains information on sites of interest ("chr", "sites", "st.posi", "en.posi", "pval.wave", "pval.ms", "qval.wave", "qval.ms", "logLR.wave", "logLR.ms", "logLR.wave.null", "logLR.ms.null")
## out.path : path to directory where figures will be saved
## wave.out.path : path to directory which contains results from wavelet analysis
## file.name : output figure file name
## siteSize : site size
## treatment : treatment name
## null : indicate whether effect size from treatment vs control (null = FALSE) or from control vs control (null = TRUE)
## strand : 'both', 'plus', 'minus'; add two strands, use + strand, or use - strand
## sig.level : +/- sig.level * standard deviation
## wave.effect : indicate whether effect size from wavelet is plotted
## multiseq.effect : indicate whether effect size from multiseq is plotted
## deseq.effect : indicate whether effect size from DESeq is plotted
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



info.path = '/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/tmp/Copper.2048.both.msOnly.info'
out.path = '/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/code_fig/JSM_2014/'
wave.out.path = '/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/wave/'
file.name= 'msOnly.1'
siteSize=2048
treatment='Copper'
null = FALSE
strand='both'
#strand='plus'
#strand='minus'
sig.level = 2
wave.effect=TRUE
multiseq.effect=TRUE
deseq.effect=FALSE




#args = (commandArgs(TRUE))
#eval(parse(text=args[[1]]))
#eval(parse(text=args[[2]]))
#eval(parse(text=args[[3]]))
#eval(parse(text=args[[4]]))
#eval(parse(text=args[[5]]))
#eval(parse(text=args[[6]]))
#eval(parse(text=args[[7]]))
#eval(parse(text=args[[8]]))
#eval(parse(text=args[[9]]))
#eval(parse(text=args[[10]]))
#eval(parse(text=args[[11]]))
#eval(parse(text=args[[12]]))

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



## read TF
all_bed = read.table(gzfile('/mnt/gluster/data/external_private_supp/roger_atacseq/dnasefootprints/OpenChromDnaseGm19239.6c.bed.gz'))
#TssAnno <- read.table('../data/Ensembl2.txt',header=F,as.is=T)
TssAnno = read.table(gzfile('/mnt/gluster/data/external_private_supp/roger_atacseq/dnasefootprints/Copper.TSS.DiffExpressed.FDR10.bed.gz'))
#/mnt/gluster/data/external_private_supp/roger_atacseq/dnasefootprints
#OpenChromDnaseGm19239.6c.bed.gz
#factorNames2.txt


## set up working directory and open figure file
setwd(out.path)
numfig = wave.effect + multiseq.effect + deseq.effect + 1
if(numfig <= 2){
    pdf(paste0(out.path, file.name, ".effect.pdf"), width=10, height=5)
}else{
    pdf(paste0(out.path, file.name, ".effect.pdf"), width=10, height=5)
}    
nf <- layout(matrix(1:numfig,numfig,1,byrow=TRUE),TRUE)


#############################
# read all information
#############################
dat.info = read.table(file=info.path, header = TRUE, as.is=TRUE)


numSites = dim(dat.info)[1]

#for(ss in 1:numSites){

ss = 1

## read location information
chr = dat.info$chr[ss]
st.posi = dat.info$st.posi[ss]
en.posi = dat.info$en.posi[ss]


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


####################
# plot raw data
####################

 
## get phenotype
xmin = st.posi
xmax = en.posi

phe.D = phenoD/library.read.depth
trt.pheno = apply(phe.D[1:(numSam/2),], 2, mean)
ctr.pheno = apply(phe.D[(numSam/2+1):numSam,], 2, mean)
trt.RC = sum(phenoD[1:(numSam/2),])
ctr.RC = sum(phenoD[(numSam/2+1):numSam,])

## smooth raw phenotype
## denoise them using wavethresh
library(wavethresh)

sig.all = c(trt.pheno, ctr.pheno)
sig.all.smooth = BAYES.THR(sig.all)
sig.all.smooth[sig.all.smooth < 0] = 0
trt.pheno.smooth = sig.all.smooth[1:siteSize]
ctr.pheno.smooth = sig.all.smooth[((1:siteSize) + siteSize)]

ymin = 0
#ymaxT = max(trt.pheno.smooth, ctr.pheno.smooth, trt.pheno, ctr.pheno)*(1+ 0.05)
ymaxT = max(trt.pheno, ctr.pheno)*(1+ 0.05)

xval = xmin:xmax
ymax = ymaxT*10^6


## get pcr information
xval_mapp = NULL
if(length(unlist(pcr.posi)) > 0){ 
    xval_mapp = xval[unlist(pcr.posi)]
}


## Make a raw phenotype figure
if(!null){
    raw.title = paste0("chr", chr, ":", st.posi, "-", en.posi, ", treatment(red):", trt.RC, " control(blue):", ctr.RC)
}else{
    raw.title = paste0("chr", chr, ":", st.posi, "-", en.posi, ", control(red):", trt.RC, " control(blue):", ctr.RC)
}    
par(mar = c(1,4,1,2))
plot(1,1,type="n", xlab = "position", ylab = "DNaseI cut rate per million reads",ylim=c(ymin, ymax),xlim=c(xmin, xmax),main =raw.title, axes=FALSE)
axis(1)
#if(!is.null(xval_mapp)){
#    axis(1, at=xval_mapp, col="green")
#}
axis(2)
box()


### Transcription factor
sel.sites = all_bed[all_bed[,1] == paste("chr", chr, sep="") & all_bed[,2] < (xmax+1) & all_bed[,3] > xmin, ]
offset = -0.0025
if(dim(sel.sites)[1] > 0){
for(k in 1:dim(sel.sites)[1]){
offset = -offset
text(x=(sel.sites[k,2] + sel.sites[k,3])/2, y=(ymax -abs(offset) - offset), strsplit(as.character(sel.sites[k,4]), split="=")[[1]][2])
rect(sel.sites[k,2], 0, sel.sites[k,3], ymax + 1, col=rgb(0,1,0,0.3), border='NA')
}
}

points(xval, ctr.pheno*10^6, col = rgb(0,0,1,alpha=0.7), type="l")
points(xval, trt.pheno*10^6, col = rgb(1,0,0,alpha=0.7), type="l")





#points(xval, trt.pheno*10^6, col = "orange", type="l")
#points(xval, ctr.pheno*10^6, col = "skyblue", type="l")
#points(xval, trt.pheno.smooth*10^6, col = "red", type="l")
#points(xval, ctr.pheno.smooth*10^6, col = "blue", type="l")


#GETS AND PLOTS ANY TSSs IN THE REGION
TSS <- TssAnno[(as.character(TssAnno[,1]) == paste("chr", chr, sep="")) & (TssAnno[,2] > xmin) & (TssAnno[,2] < (xmax+1)),]
if(dim(TSS)[1] > 0) {
for(k in 1:dim(TSS)[1]){
mtext('*', side=1, at=TSS[k,2], col='purple', cex=1.5, padj=1)
}
}







########################
# multiseq effect size
########################

if(multiseq.effect){

    ## title
    if(!null){
        title = paste0("multiseq [+/-", sig.level, "] -log10(pval): ", round(-log(dat.info$pval.ms[ss],10),2), " -log10(qval): ", round(-log(dat.info$qval.ms[ss],10),2), " logLR: ", round(dat.info$logLR.ms[ss],2))
    }else{
        title = paste0("multiseq [+/-", sig.level, "] logLR: ", round(dat.info$logLR.ms.null[ss],2))
    }
        
    ## get effect size
    genoD = g
    res = multiseq(x = phenoD, g = genoD, read.depth = library.read.depth)

    effect.mean = -res$effect.mean
    effect.sd = sqrt(res$effect.var)

    effect.low = effect.mean - sig.level*effect.sd
    effect.high= effect.mean + sig.level*effect.sd

    ymax = max(effect.high) + 10^(-7)
    ymin = min(effect.low) - 10^(-7)
    
    wh.low = which(effect.low > 0)
    wh.high = which(effect.high < 0)
    high_wh = sort(unique(union(wh.low, wh.high)))
    col_posi = xval[high_wh]

    par(mar = c(1,4,4,2))
    plot(1,1,type="n", xlab = "position", ylab = "Effect size",ylim=c(ymin, ymax),xlim=c(xmin, xmax),main =title, axes=FALSE)
    axis(2)
    if(length(col_posi) > 0){
        for(j in 1:length(col_posi)){
            polygon(c(col_posi[j]-0.5, col_posi[j]-0.5, col_posi[j]+0.5, col_posi[j]+0.5), c(ymin-2, ymax+2, ymax+2, ymin-2), col ="pink", border = NA)
        }
    }
    abline(h = 0, col = "red")
    points(xval, effect.mean, col = "blue", type="l")
    points(xval, effect.high, col = "skyblue", type="l")
    points(xval, effect.low, col = "skyblue", type="l")
    box()

}





########################
# wavelet effect size
########################

if(wave.effect){

    ## title
    if(!null){
        title = paste0("wavelet [+/-", sig.level, "] -log10(pval): ", round(-log(dat.info$pval.wave[ss],10),2), " -log10(qval): ", round(-log(dat.info$qval.wave[ss],10),2), " logLR: ", round(dat.info$logLR.wave[ss],2))
    }else{
        title = paste0("wavelet [+/-", sig.level, "] logLR: ", round(dat.info$logLR.wave.null[ss],2))
    }
    
    ## get effect size
    Wmat = read.table(paste0(WaveQTL.repodir, "data/DWT/Wmat_", siteSize), as.is = TRUE)
    W2mat = Wmat*Wmat

    if(!null){
        path.wave = paste0(wave.out.path, treatment, ".", siteSize, ".", strand, ".alt.run/output/", treatment, ".", siteSize, ".", strand, ".alt.", chr, ".", dat.info$sites[ss], ".fph.")
    }else{
        path.wave = paste0(wave.out.path, treatment, ".", siteSize, ".", strand, ".null.run/output/", treatment, ".", siteSize, ".", strand, ".null.", chr, ".", dat.info$sites[ss], ".fph.")
    }

    effect.mean.w = as.numeric(read.table(paste0(path.wave, "mean.txt", sep=""))[-1])
    effect.mean = -matrix(data=effect.mean.w, nr = 1, nc = siteSize)%*%as.matrix(Wmat)

    effect.var.w = as.numeric(read.table(paste0(path.wave, "var.txt", sep=""))[-1])
    effect.sd = sqrt(matrix(data=effect.var.w, nr = 1, nc = siteSize)%*%as.matrix(W2mat))

    effect.low = effect.mean - sig.level*effect.sd
    effect.high= effect.mean + sig.level*effect.sd

    ymax = max(effect.high) + 10^(-7)
    ymin = min(effect.low) - 10^(-7)
    
    wh.low = which(effect.low > 0)
    wh.high = which(effect.high < 0)
    high_wh = sort(unique(union(wh.low, wh.high)))
    col_posi = xval[high_wh]

    par(mar = c(1,4,2,2))
    plot(1,1,type="n", xlab = "position", ylab = "Effect size",ylim=c(ymin, ymax),xlim=c(xmin, xmax),main =title, axes=FALSE)
    axis(2)
    if(length(col_posi) > 0){
        for(j in 1:length(col_posi)){
            polygon(c(col_posi[j]-0.5, col_posi[j]-0.5, col_posi[j]+0.5, col_posi[j]+0.5), c(ymin-2, ymax+2, ymax+2, ymin-2), col ="pink", border = NA)
        }
    }
    abline(h = 0, col = "red")
    points(xval, effect.mean, col = "blue", type="l")
    points(xval, effect.high, col = "skyblue", type="l")
    points(xval, effect.low, col = "skyblue", type="l")
    box()

}


dev.off()


 

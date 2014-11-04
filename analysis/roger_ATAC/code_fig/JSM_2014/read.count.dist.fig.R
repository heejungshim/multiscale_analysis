#!/usr/bin/env Rscript

## Aim : This file contains Rscript to make figure of read count distribution
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+



setwd("/mnt/lustre/home/shim/multiscale_analysis")
multiscale.analysis.repodir <- scan(".multiscale_analysis.repodir.txt", what=character())


##############
# msOnly
##############

info.path = '/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/tmp/Copper.2048.both.msOnly.info'
siteSize=2048
treatment='Copper'
null = FALSE
strand='both'
deseq.dat.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/"
numSam = 6
numS = siteSize %/% 300

## read all information

dat.info = read.table(file=info.path, header = TRUE, as.is=TRUE)
numSites = dim(dat.info)[1]

rc.dat = matrix(data=NA, nr = numSites, nc = numSam)

for(chr in 1:22){
    wh.here = which(dat.info$chr == chr)
    if(length(wh.here) > 0){
        ## read count data
        deseq.dat = read.table(paste0(deseq.dat.path, treatment, ".", siteSize, ".", strand, ".300.alt.run/data.each.site.", chr, ".txt"), as.is = TRUE)
        rc.dat[wh.here,] = as.matrix(deseq.dat[dat.info$sites[wh.here],])
    }
}



msOnly.rc = apply(rc.dat,1,sum)




####################
# both ms and wave
####################

info.path = '/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/tmp/Copper.2048.both.both.info'
siteSize=2048
treatment='Copper'
null = FALSE
strand='both'
deseq.dat.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/"
numSam = 6
numS = siteSize %/% 300

# read all information
dat.info = read.table(file=info.path, header = TRUE, as.is=TRUE)
numSites = dim(dat.info)[1]

rc.dat = matrix(data=NA, nr = numSites, nc = numSam)

for(chr in 1:22){
    wh.here = which(dat.info$chr == chr)
    if(length(wh.here) > 0){
        ## read count data
        deseq.dat = read.table(paste0(deseq.dat.path, treatment, ".", siteSize, ".", strand, ".300.alt.run/data.each.site.", chr, ".txt"), as.is = TRUE)
        rc.dat[wh.here,] = as.matrix(deseq.dat[dat.info$sites[wh.here],])
    }
}


both.rc = apply(rc.dat,1,sum)


#########################
# all sites
#########################



siteSize=2048
treatment='Copper'
null = FALSE
strand='both'
deseq.dat.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/"
numSam = 6
numS = siteSize %/% 300

## read total number of sites
path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/", treatment, ".", siteSize, ".numSites.txt")
numSites.list = scan(path)

rc.dat = vector("list", 22)
for(chr in 1:22){
    deseq.dat = read.table(paste0(deseq.dat.path, treatment, ".", siteSize, ".", strand, ".300.alt.run/data.each.site.", chr, ".txt"), as.is = TRUE)
    rc.dat[[chr]] = apply(deseq.dat,1, sum)
}

temp = unlist(rc.dat)
all.rc = temp[temp > 0]



median(msOnly.rc)
median(both.rc)
median(all.rc)

#> median(msOnly.rc)
#[1] 409
#> median(both.rc)
#[1] 1796
#> median(all.rc)
#[1] 239






rc = c(msOnly.rc, both.rc, all.rc)
cate = c(rep("only multiseq", length(msOnly.rc)), rep("multiseq and WaveQTL", length(both.rc)), rep("all sites", length(all.rc)))
wh = which(rc > 7000)
setwd("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/code_fig/JSM_2014/")
pdf("distReadCount.pdf")
boxplot(rc[-wh]~cate[-wh], main = "", ylab = "read count")
dev.off()


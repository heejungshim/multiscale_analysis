## `see.read.count.dist.R' contains script to see distribution of read count on two extreme cases...
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
cate = c(rep("only multiseq", length(msOnly.rc)), rep("multiseq and wave", length(both.rc)), rep("all", length(all.rc)))
wh = which(rc > 300000)
pdf("distReadCount.new.pdf")
boxplot(rc[-wh]~cate[-wh], main = "", ylab = "read count")
dev.off()


#####################
# FDR = 0.01 in ms
#####################


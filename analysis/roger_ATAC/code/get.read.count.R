## `get.read.count.R' contains script to get read count for each site
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




setwd("/mnt/lustre/home/shim/multiscale_analysis")
multiscale.analysis.repodir <- scan(".multiscale_analysis.repodir.txt", what=character())


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



for(chr in 1:22){
    numSites = numSites.list[chr]
    ## read count data
    deseq.dat = read.table(paste0(deseq.dat.path, treatment, ".", siteSize, ".", strand, ".300.alt.run/data.", chr, ".txt"), as.is = TRUE)
    st = (1:numSites)*numS - (numS -1)
    en = (1:numSites)*numS
    res = sapply(1:numSites, function(x, st.in, en.in, deseq.dat.in){return(as.numeric(apply(deseq.dat.in[(st[x]):(en[x]),],2,sum)))}, st.in = st, en.in=en, deseq.dat.in = deseq.dat)

    path.out = paste0(deseq.dat.path, treatment, ".", siteSize, ".", strand, ".300.alt.run/data.each.site.", chr, ".txt")
    write.table(t(res), file=path.out, quote=FALSE, row.names=FALSE, col.names = FALSE)
}
    

 



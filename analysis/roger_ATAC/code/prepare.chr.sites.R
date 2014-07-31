#!/usr/bin/env Rscript

## Aim : This file contains Rscripts to prepare chromosome and sites corresponding to list of test statistics.
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+

siteSize=2048
treatment='Copper'

## read number of sites for each chromosome
path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/", treatment, ".", siteSize, ".numSites.txt")
numSites.list = scan(path)

chr.list = sites.list = vector("list", 22)
for(chr in 1:22){
    chr.list[[chr]] = rep(chr, numSites.list[chr])
    sites.list[[chr]] = 1:numSites.list[chr]
}

out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/info/", treatment, ".", siteSize, ".chr.sites.Robj")

save("chr.list", "sites.list", file = out.path)






siteSize=1024
treatment='Copper'

## read number of sites for each chromosome
path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/", treatment, ".", siteSize, ".numSites.txt")
numSites.list = scan(path)

chr.list = sites.list = vector("list", 22)
for(chr in 1:22){
    chr.list[[chr]] = rep(chr, numSites.list[chr])
    sites.list[[chr]] = 1:numSites.list[chr]
}

out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/info/", treatment, ".", siteSize, ".chr.sites.Robj")

save("chr.list", "sites.list", file = out.path)






siteSize=2048
treatment='Selenium'

## read number of sites for each chromosome
path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/", treatment, ".", siteSize, ".numSites.txt")
numSites.list = scan(path)

chr.list = sites.list = vector("list", 22)
for(chr in 1:22){
    chr.list[[chr]] = rep(chr, numSites.list[chr])
    sites.list[[chr]] = 1:numSites.list[chr]
}

out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/info/", treatment, ".", siteSize, ".chr.sites.Robj")

save("chr.list", "sites.list", file = out.path)





siteSize=1024
treatment='Selenium'

## read number of sites for each chromosome
path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/", treatment, ".", siteSize, ".numSites.txt")
numSites.list = scan(path)

chr.list = sites.list = vector("list", 22)
for(chr in 1:22){
    chr.list[[chr]] = rep(chr, numSites.list[chr])
    sites.list[[chr]] = 1:numSites.list[chr]
}

out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/info/", treatment, ".", siteSize, ".chr.sites.Robj")

save("chr.list", "sites.list", file = out.path)






siteSize=2048
treatment='Retinoic'

## read number of sites for each chromosome
path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/", treatment, ".", siteSize, ".numSites.txt")
numSites.list = scan(path)

chr.list = sites.list = vector("list", 22)
for(chr in 1:22){
    chr.list[[chr]] = rep(chr, numSites.list[chr])
    sites.list[[chr]] = 1:numSites.list[chr]
}

out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/info/", treatment, ".", siteSize, ".chr.sites.Robj")

save("chr.list", "sites.list", file = out.path)







siteSize=1024
treatment='Retinoic'

## read number of sites for each chromosome
path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/", treatment, ".", siteSize, ".numSites.txt")
numSites.list = scan(path)

chr.list = sites.list = vector("list", 22)
for(chr in 1:22){
    chr.list[[chr]] = rep(chr, numSites.list[chr])
    sites.list[[chr]] = 1:numSites.list[chr]
}

out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/info/", treatment, ".", siteSize, ".chr.sites.Robj")

save("chr.list", "sites.list", file = out.path)



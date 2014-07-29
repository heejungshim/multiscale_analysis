#!/usr/bin/env Rscript

## Aim : This file contains Rscript to collect logLR from different chromosome 
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+



wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/wave/'
siteSize=2048
treatment='Copper'
null=FALSE
strand='both'

    
## directory name 
if(!null){
    com.dir.name=paste0(treatment,".", siteSize, ".", strand, ".alt")
}else{
    com.dir.name=paste0(treatment,".", siteSize, ".", strand, ".null")
}


path = paste0(wd.path, com.dir.name, ".run/sum/", com.dir.name, ".")
done_list = vector("list", 22)
logLR_list = vector("list", 22)

for(chr in 1:22){
    load(paste0(path, chr, ".Robj"))
    done_list[[chr]] = done
    logLR_list[[chr]] = logLR
}


save("logLR_list", "done_list", file =paste0(path, "Robj"))




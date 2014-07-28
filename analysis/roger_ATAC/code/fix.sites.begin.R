#!/usr/bin/env Rscript

## Aim : This file contains Rscripts to fix the first two sites if they are close to the start of chromosome
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+



#Treatment = "Selenium"
Treatment = "Retinoic"
#Treatment = "Copper"

#size = 2048
size = 1024

## path to location information 
path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/"

## read chromosome length
path.chr.len = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/info/"
chr.len.list = scan(file = paste0(path.chr.len, "chr.len.txt"))

case.two = NA

problem.site = vector("list", 22)
for(chr in 1:22){
    res = read.table(file = paste0(path, Treatment, ".", size, ".chr", chr, ".locus"),header= TRUE) 
    chr.len = chr.len.list[chr]
    wh = which(res$st.posi  < 5)
    problem.site[[chr]] = wh
    
    if(length(wh) > 0){
        if(length(wh) > 1){
            case.two = TRUE
        }

        res$st.posi[wh] = 5
        res$en.posi[wh] = 5 + size  - 1
        res = data.frame(chr = res$chr, st.posi = res$st.posi, en.posi = res$en.posi)
        write.table(res, file = paste0(path, Treatment, ".", size, ".chr", chr, ".locus"), col.names = TRUE, row.names= FALSE, quote= FALSE) 
    }
    
}

case.two
problem.site






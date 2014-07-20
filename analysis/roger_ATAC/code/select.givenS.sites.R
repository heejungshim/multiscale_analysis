#!/usr/bin/env Rscript

## Aim : This file contains Rscripts to combine two adjacent top 5% 300bp windows that are within a given size. 
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+



#Treatment = "Selenium"
Treatment = "Retinoic"
#Treatment = "Copper"

#size = 2048
size = 1024


for(chr in 1:22){


# read top 5% of 300bp region
path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/info/"
top5.sites = read.table(file = paste0(path, Treatment, ".05.txt"), as.is = TRUE)


# focus on given chromosome
wh = which(top5.sites[,1] == paste0("chr", chr))
sel.top5 = top5.sites[wh,]

# get chromosome length
path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/info/"
chr.len = scan(file = paste0(path, "chr.len.txt"))[chr]

# combine adjacent 300bp windows
len = dim(sel.top5)[1]
st = 1
center.posi = NULL
center.posi.ix = 1
while(st <= len){
    st.posi = sel.top5[st,2]
    max.posi = st.posi + size - 1
    wh = which(sel.top5[,3] <= max.posi)
    en = wh[length(wh)]
    en.posi = sel.top5[en,3]
    center.posi[center.posi.ix] = (en.posi + st.posi)/2
    center.posi.ix = center.posi.ix + 1 
    st = en + 1
}

if(center.posi[1] - size/2 < 1){
    center.posi[1] = size/2 + 2
}

if(center.posi[length(center.posi)] + size/2 > chr.len){
    center.posi[length(center.posi)] = chr.len - size/2 - 2
}


path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/"
res = data.frame(chr = rep(paste0("chr",chr), length(center.posi)), st.posi = center.posi - size/2 + 1, en.posi = center.posi + size/2)
write.table(res, file = paste0(path, Treatment, ".", size, ".chr", chr, ".locus"), col.names = TRUE, row.names= FALSE, quote= FALSE) 


}

## check if I did correctly!
#len = dim(sel.top5)[1]
#err = NULL
#err.ix = 1
#for(i in 1:len){
#    i = 1
#    wh = which((sel.top5[i,2] >= res$st.posi) & (sel.top5[i,3] <= res$en.posi))
#    if(length(wh) == 0){
#        err[err.ix] = i
#        err.ix = err.ix + 1
#    }
#}

#err
# NULL

#!/usr/bin/env Rscript

## Aim : This file contains Rscripts to get total number of sites (1024 or 2048 bps) for each chromosome.
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+



#Treatment = "Selenium"
Treatment = "Retinoic"
#Treatment = "Copper"

#size = 2048
size = 1024

path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/"

total.num = rep(NA, 22)
for(chr in 1:22){
    res = read.table(file = paste0(path, Treatment, ".", size, ".chr", chr, ".locus"),header= TRUE)
    total.num[chr] = dim(res)[1]
}



output.path = paste0(path, Treatment, ".", size, ".numSites.txt")
cat(total.num, file = output.path)
sum(total.num)



# "Copper" 1024 :  258026
# "Copper" 2048 :  197969
# "Selenium" 1024 : 260742
# "Selenium" 2048 : 200352
# "Retinoic" 1024 : 259761
# "Retinoic" 2048 : 197506



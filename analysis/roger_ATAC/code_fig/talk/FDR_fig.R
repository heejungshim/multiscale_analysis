#!/usr/bin/env Rscript

## Aim : This file contains Rscript to make figure for FDR curves
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+



## copy from see.statistic.3method.both.R


siteSize=2048
treatment='Copper'
strand='both'


all.name = paste0(treatment,".", siteSize, ".", strand)



## read p-value

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/deseq.", all.name, ".Robj"))

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/wave.new.", all.name, ".Robj"))

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/ms.new.", all.name, ".Robj"))




#####################
### Storey FDR 
#####################

#install.packages("qvalue")
library("qvalue")

qval.wave = qvalue(pval.wave.new)
qval.ms = qvalue(pval.ms.new)
qval.deseq.100 = qvalue(pval.deseq.100)
qval.deseq.300 = qvalue(pval.deseq.300)

qval.wave$pi0
qval.ms$pi0
qval.deseq.100$pi0
qval.deseq.300$pi0

#> qval.wave$pi0
#[1] 0.7818706
#> qval.ms$pi0
#[1] 0.883166
#> qval.deseq.100$pi0
#[1] 0.9467088
#> qval.deseq.300$pi0
#[1] 0.9528259


alpha.list = seq(0, 0.2, by=0.01)
length(alpha.list)
# 20
num.wave = num.ms = num.deseq.100 = num.deseq.300 = rep(NA, length(alpha.list))
for(i in 1:length(alpha.list)){
    num.wave[i] = sum(qval.wave$qvalues < alpha.list[i])
    num.ms[i] = sum(qval.ms$qvalues < alpha.list[i])
    num.deseq.100[i] = sum(qval.deseq.100$qvalues < alpha.list[i])
    num.deseq.300[i] = sum(qval.deseq.300$qvalues < alpha.list[i])
}



setwd("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/code/")
pdf("discovery_FDR_JSM2014.pdf")

ymax = max(num.wave, num.ms, num.deseq.100, num.deseq.300)
ymin = 0
plot(alpha.list, num.ms, ylim=c(ymin,ymax), col="red", type = "l", ylab = "number of significant tests", xlab = "FDR", main = "", lty=1, cex=1)
points(alpha.list, num.wave, ylim=c(ymin,ymax), col="blue", type="l",lty=1, cex=1)
#points(alpha.list, num.deseq.100, ylim=c(ymin,ymax), col="darkgreen", type="l",lty=1, cex=1)
#points(alpha.list, num.deseq.300, ylim=c(ymin,ymax), col="darkgreen", type="l", lty="dashed",lty=1, cex=1)

legend(0,ymax, c("multiseq", "WaveQTL"), col = c("red",  "blue"), lty = c(1,1), cex = 1, text.col = "black",merge = FALSE, bg = "white")

dev.off()


###########################
# check intersection
###########################


alpha = 0.2

sum(qval.wave$qvalues < alpha)
sum(qval.ms$qvalues < alpha)
sum((qval.wave$qvalues < alpha) & (qval.ms$qvalues < alpha))

sum((qval.wave$qvalues < alpha) & (qval.ms$qvalues < alpha))/sum(qval.ms$qvalues < alpha)
sum((qval.wave$qvalues < alpha) & (qval.ms$qvalues < alpha))/sum(qval.wave$qvalues < alpha)


#> alpha = 0.15
#> 
#> sum(qval.wave$qvalues < alpha)
#[1] num_wave[1] = 236
#> sum(qval.ms$qvalues < alpha)
#[1] num_multiseq[1] = 2293
#> sum((qval.wave$qvalues < alpha) & (qval.ms$qvalues < alpha))
#[1] num_both[1] =  165
#> 
#> sum((qval.wave$qvalues < alpha) & (qval.ms$qvalues < alpha))/sum(qval.ms$qvalues < alpha)
#[1] 0.07195813
#> sum((qval.wave$qvalues < alpha) & (qval.ms$qvalues < alpha))/sum(qval.wave$qvalues < alpha)
#[1] 0.6991525


#> alpha = 0.2
#> 
#> sum(qval.wave$qvalues < alpha)
#[1] 377
#> sum(qval.ms$qvalues < alpha)
#[1] 3477
#> sum((qval.wave$qvalues < alpha) & (qval.ms$qvalues < alpha))
#[1] 301
#> 
#> sum((qval.wave$qvalues < alpha) & (qval.ms$qvalues < alpha))/sum(qval.ms$qvalues < alpha)
#[1] 0.08656888
#> sum((qval.wave$qvalues < alpha) & (qval.ms$qvalues < alpha))/sum(qval.wave$qvalues < alpha)
#[1] 0.7984085

#!/usr/bin/env Rscript

## Aim : make FDR curves from WaveQTL, multiseq, and DESeq2 with different bin sizes. I modified a script here (cp ~/multiscale_analysis/analysis/roger_ATAC/code/see.FDR.R).
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+



siteSize=2048
treatment='Copper'
strand='both'
window.size=300
numSam = 6
all.name = paste0(treatment,".", siteSize, ".", strand)


out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.pval.discrete.Robj")

load(out.path)
## pval.deseq.3.0
## pval.deseq.3.60
## pval.deseq.3.30
## pval.deseq.3.20
## pval.deseq.3.10



out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".alt.run/output/DESeq2.pval.discrete.",600,".Robj")

load(out.path)


#####################
### Storey FDR 
#####################

#install.packages("qvalue")
library("qvalue")

input.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/pval.ms.wave.", all.name, ".Robj")
load(input.path)
## pval.ms pval.wave


input.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".alt.run/output/DESeq2.pval.discrete.", 10, ".Robj")
load(input.path)

input.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".alt.run/output/DESeq2.pval.discrete.Robj")
load(input.path)

#######################
## cut.thresh = 0
## 100, 300, 600, 2048
#######################


## "pval.deseq.100.0"  "pval.deseq.100.10" "pval.deseq.3.0"    "pval.deseq.3.10"  "pval.deseq.full.0" "pval.ms" "pval.wave"  

sum(is.na(pval.deseq.600.0))
sum(is.na(pval.deseq.100.0))
sum(is.na(pval.deseq.100.10))
sum(is.na(pval.deseq.3.0))
sum(is.na(pval.deseq.3.10))
sum(is.na(pval.deseq.full.0))
sum(is.na(pval.ms))
sum(is.na(pval.wave))

## [1] 11565
## [1] 11565
## [1] 16805
## [1] 11565
## [1] 13014
## [1] 11565
## [1] 11570
## [1] 11570

length(pval.ms)
## [1] 197969


del.ix = NULL
del.ix = union(del.ix, which(is.na(pval.deseq.100.0)==TRUE))
del.ix = union(del.ix, which(is.na(pval.deseq.3.0)==TRUE))
del.ix = union(del.ix, which(is.na(pval.deseq.600.0)==TRUE))
del.ix = union(del.ix, which(is.na(pval.deseq.full.0)==TRUE))
del.ix = union(del.ix, which(is.na(pval.ms)==TRUE))
del.ix = union(del.ix, which(is.na(pval.wave)==TRUE))

length(del.ix)
## 11570

qval.wave = qvalue(pval.wave[-del.ix])
qval.ms = qvalue(pval.ms[-del.ix])
qval.deseq.full.0 = qvalue(pval.deseq.full.0[-del.ix])
qval.deseq.3.0 = qvalue(pval.deseq.3.0[-del.ix])
qval.deseq.100.0 = qvalue(pval.deseq.100.0[-del.ix])
qval.deseq.600.0 = qvalue(pval.deseq.600.0[-del.ix])



qval.wave$pi0
qval.ms$pi0
qval.deseq.full.0$pi0
qval.deseq.3.0$pi0
qval.deseq.100.0$pi0
qval.deseq.600.0$pi0

## [1] 0.7818706
## [1] 0.883166
## [1] 0.9649561
## [1] 0.9448775
## [1] 1
## [1] 0.9551462

alpha.list = seq(0, 0.2, by=0.01)
length(alpha.list)
# 20
num.wave = num.ms = num.deseq.full.0 = num.deseq.3.0 = num.deseq.100.0 = num.deseq.600.0 = rep(NA, length(alpha.list))
for(i in 1:length(alpha.list)){
    num.wave[i] = sum(qval.wave$qvalues < alpha.list[i])
    num.ms[i] = sum(qval.ms$qvalues < alpha.list[i])
    num.deseq.full.0[i] = sum(qval.deseq.full.0$qvalues < alpha.list[i])
    num.deseq.3.0[i] = sum(qval.deseq.3.0$qvalues < alpha.list[i])
    num.deseq.100.0[i] = sum(qval.deseq.100.0$qvalues < alpha.list[i])
    num.deseq.600.0[i] = sum(qval.deseq.600.0$qvalues < alpha.list[i])
}




pdf("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/code_fig/RECOMB_2014/fig/FDR.pdf")


ymax = max(num.wave, num.ms, num.deseq.full.0, num.deseq.3.0, num.deseq.100.0, num.deseq.600.0)
ymin = 0

plot(alpha.list, num.ms, ylim=c(ymin,ymax), col="red", type = "l", main="", xlab = "FDR", ylab="number of significant tests")
points(alpha.list, num.wave, ylim=c(ymin,ymax), col="blue", type="l")
points(alpha.list, num.deseq.full.0, ylim=c(ymin,ymax), col="darkgreen", type="l")
points(alpha.list, num.deseq.3.0, ylim=c(ymin,ymax), col="orange", type="l")
points(alpha.list, num.deseq.100.0, ylim=c(ymin,ymax), col="purple", type="l")
points(alpha.list, num.deseq.600.0, ylim=c(ymin,ymax), col="brown", type="l")
legend(0.01, ymax*0.99, c("multiseq", "WaveQTL", "DESeq2 2048bp", "DESeq2 600bp", "DESeq2 300bp", "DESeq2 100bp") , col = c("red", "blue", "darkgreen", "brown", "orange", "purple"), lty = c(1,1,1, 1, 1, 1), text.col = "black",merge = FALSE, bg = "white")


dev.off()




#!/usr/bin/env Rscript

## Aim : check if results from different approaches overlap by using q-values. It wasn't informative. So I checked number of overlapping detection. 
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


output.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/fig/CompareMethods/overlap.pdf"

pdf(output.path)


smoothScatter(-log(qval.wave$qvalues,10), -log(qval.ms$qvalues,10), xlab = "WaveQTL", ylab = "multiseq", main = "-log10(qvalues)")
points(c(0,6), c(0,6), col="red", type="l")

smoothScatter(-log(qval.deseq.full.0$qvalues,10), -log(qval.ms$qvalues,10), xlab = "DESeq2 2048", ylab = "multiseq", main = "-log10(qvalues)")
points(c(0,6), c(0,6), col="red", type="l")

smoothScatter(-log(qval.deseq.3.0$qvalues,10), -log(qval.ms$qvalues,10), xlab = "DESeq2 300", ylab = "multiseq", main = "-log10(qvalues)")
points(c(0,6), c(0,6), col="red", type="l")

smoothScatter(-log(qval.deseq.3.0$qvalues,10), -log(qval.deseq.full.0$qvalues,10), xlab = "DESeq2 300", ylab = "DESeq2 2048", main = "-log10(qvalues)")
points(c(0,6), c(0,6), col="red", type="l")

smoothScatter(-log(qval.deseq.600.0$qvalues,10), -log(qval.deseq.full.0$qvalues,10), xlab = "DESeq2 600", ylab = "DESeq2 2048", main = "-log10(qvalues)")
points(c(0,6), c(0,6), col="red", type="l")

smoothScatter(-log(qval.deseq.3.0$qvalues,10), -log(qval.deseq.600.0$qvalues,10), xlab = "DESeq2 300", ylab = "DESeq2 600", main = "-log10(qvalues)")
points(c(0,6), c(0,6), col="red", type="l")



##plot(-log(qval.wave$qvalues,10), -log(qval.ms$qvalues,10), xlab = "WaveQTL", ylab = "multiseq", main = "-log10(qvalues)")
##plot(-log(qval.deseq.full.0$qvalues,10), -log(qval.ms$qvalues,10), xlab = "DESeq2 2048", ylab = "multiseq", main = "-log10(qvalues)")
##plot(-log(qval.deseq.3.0$qvalues,10), -log(qval.ms$qvalues,10), xlab = "DESeq2 300", ylab = "multiseq", main = "-log10(qvalues)")
##plot(-log(qval.deseq.3.0$qvalues,10), -log(qval.deseq.full.0$qvalues,10), xlab = "DESeq2 300", ylab = "DESeq2 2048", main = "-log10(qvalues)")
##plot(-log(qval.deseq.600.0$qvalues,10), -log(qval.deseq.full.0$qvalues,10), xlab = "DESeq2 600", ylab = "DESeq2 2048", main = "-log10(qvalues)")
##plot(-log(qval.deseq.3.0$qvalues,10), -log(qval.deseq.600.0$qvalues,10), xlab = "DESeq2 300", ylab = "DESeq2 600", main = "-log10(qvalues)")


dev.off()


##############################################
##
## Let's look at percentage of overlapping at different FDR (0.05, 0.1, 0.2)
##
##############################################

cut = c(0.05, 0.1, 0.15)

q1 = qval.wave$qvalues
q2 = qval.ms$qvalues

num1 = num2 = num12  = rep(NA, length(cut))

for(i in 1:length(cut)){
  num1[i] = sum(q1 < cut[i])
  num2[i] = sum(q2 < cut[i])
  num12[i] = sum((q1 < cut[i]) & (q2 < cut[i]))
}

num1
num2
num12
num12/num1
num12/num2


## [1]   0  16 236
## [1]  819 1428 2293
## [1]   0  11 165
## [1]       NaN 0.6875000 0.6991525
## [1] 0.000000000 0.007703081 0.071958133


q1 = qval.deseq.full.0$qvalues
q2 = qval.ms$qvalues

num1 = num2 = num12  = rep(NA, length(cut))

for(i in 1:length(cut)){
  num1[i] = sum(q1 < cut[i])
  num2[i] = sum(q2 < cut[i])
  num12[i] = sum((q1 < cut[i]) & (q2 < cut[i]))
}

num1
num2
num12
num12/num1
num12/num2


## [1] 322 530 683
## [1]  819 1428 2293
## [1] 190 335 467
## [1] 0.5900621 0.6320755 0.6837482
## [1] 0.2319902 0.2345938 0.2036633




q1 = qval.deseq.600.0$qvalues
q2 = qval.ms$qvalues

num1 = num2 = num12  = rep(NA, length(cut))

for(i in 1:length(cut)){
  num1[i] = sum(q1 < cut[i])
  num2[i] = sum(q2 < cut[i])
  num12[i] = sum((q1 < cut[i]) & (q2 < cut[i]))
}

num1
num2
num12
num12/num1
num12/num2

## [1]  96 260 397
## [1]  819 1428 2293
## [1]  67 180 280
## [1] 0.6979167 0.6923077 0.7052897
## [1] 0.08180708 0.12605042 0.12211077


q1 = qval.deseq.3.0$qvalues
q2 = qval.ms$qvalues

num1 = num2 = num12  = rep(NA, length(cut))

for(i in 1:length(cut)){
  num1[i] = sum(q1 < cut[i])
  num2[i] = sum(q2 < cut[i])
  num12[i] = sum((q1 < cut[i]) & (q2 < cut[i]))
}

num1
num2
num12
num12/num1
num12/num2

## [1] 286 442 729
## [1]  819 1428 2293
## [1] 143 245 404
## [1] 0.5000000 0.5542986 0.5541838
## [1] 0.1746032 0.1715686 0.1761884




q1 = qval.deseq.3.0$qvalues
q2 = qval.deseq.full.0$qvalues

num1 = num2 = num12  = rep(NA, length(cut))

for(i in 1:length(cut)){
  num1[i] = sum(q1 < cut[i])
  num2[i] = sum(q2 < cut[i])
  num12[i] = sum((q1 < cut[i]) & (q2 < cut[i]))
}

num1
num2
num12
num12/num1
num12/num2

## [1] 286 442 729
## [1] 322 530 683
## [1] 168 252 348
## [1] 0.5874126 0.5701357 0.4773663
## [1] 0.5217391 0.4754717 0.5095168




q1 = qval.deseq.600.0$qvalues
q2 = qval.deseq.full.0$qvalues

num1 = num2 = num12  = rep(NA, length(cut))

for(i in 1:length(cut)){
  num1[i] = sum(q1 < cut[i])
  num2[i] = sum(q2 < cut[i])
  num12[i] = sum((q1 < cut[i]) & (q2 < cut[i]))
}

num1
num2
num12
num12/num1
num12/num2

## [1]  96 260 397
## [1] 322 530 683
## [1]  81 184 261
## [1] 0.8437500 0.7076923 0.6574307
## [1] 0.2515528 0.3471698 0.3821376


q1 = qval.deseq.3.0$qvalues
q2 = qval.deseq.600.0$qvalues

num1 = num2 = num12  = rep(NA, length(cut))

for(i in 1:length(cut)){
  num1[i] = sum(q1 < cut[i])
  num2[i] = sum(q2 < cut[i])
  num12[i] = sum((q1 < cut[i]) & (q2 < cut[i]))
}

num1
num2
num12
num12/num1
num12/num2

## [1] 286 442 729
## [1]  96 260 397
## [1]  90 197 309
## [1] 0.3146853 0.4457014 0.4238683
## [1] 0.9375000 0.7576923 0.7783375



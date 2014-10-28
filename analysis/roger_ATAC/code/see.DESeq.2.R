#!/usr/bin/env Rscript

## Aim : Resutls from DESeq analysis are very bad. I'm not sure I used DESseq properly. So I'll try different ways to use DESeq and see if results make sens. For this purpose, I'll focus on `Copper.2048.both.300'. I'll try to exclude 300bp windows with low read counts.
## In particular, here see.DESeq.2.R I read outputs from see.DESeq.R and make figures (histogram and fdr..etc)
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+



######################################################
##
## Look at results from DESeq with different cutoff
##
######################################################

siteSize=2048
treatment='Copper'
strand='both'
window.size=300
numSam = 6
all.name = paste0(treatment,".", siteSize, ".", strand)


out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/deseq.pval.Robj")

load(out.path)
## pval.deseq.3
## pval.deseq.3.60
## pval.deseq.3.30
## pval.deseq.3.20
## pval.deseq.3.10





#################################################
## make a histogram of statistics from null and alternative and p-values for DESeq
#################################################

setwd("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/fig/DESeqdebug/")
pdf("hist.statistic.pval.DESeq.debug.pdf")

#### filter.cut = 0

filter.cut = 0
pval = pval.deseq.3

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.txt"))[,1])
 
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))

xmax = 1
xmin = 0

par(mfrow = c(3,1))
hist(deseq.alt[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : alt test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(deseq.null[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : null test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(pval[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : pvalue ", length(pval)-length(del.ix.deseq)), xlim=c(xmin, xmax), xlab ="p=value")


#### filter.cut = 10

filter.cut = 10
pval = pval.deseq.3.10

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.", filter.cut, ".txt"))[,1])
 
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))

xmax = 1
xmin = 0

par(mfrow = c(3,1))
hist(deseq.alt[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : alt test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(deseq.null[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : null test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(pval[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : pvalue ", length(pval)-length(del.ix.deseq)), xlim=c(xmin, xmax), xlab ="p=value")




#### filter.cut = 20

filter.cut = 20
pval = pval.deseq.3.20

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.", filter.cut, ".txt"))[,1])
 
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))

xmax = 1
xmin = 0

par(mfrow = c(3,1))
hist(deseq.alt[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : alt test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(deseq.null[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : null test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(pval[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : pvalue ", length(pval)-length(del.ix.deseq)), xlim=c(xmin, xmax), xlab ="p=value")



#### filter.cut = 30

filter.cut = 30
pval = pval.deseq.3.30

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.", filter.cut, ".txt"))[,1])
 
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))

xmax = 1
xmin = 0

par(mfrow = c(3,1))
hist(deseq.alt[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : alt test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(deseq.null[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : null test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(pval[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : pvalue ", length(pval)-length(del.ix.deseq)), xlim=c(xmin, xmax), xlab ="p=value")





#### filter.cut = 60

filter.cut = 60
pval = pval.deseq.3.60

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.", filter.cut, ".txt"))[,1])
 
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))

xmax = 1
xmin = 0

par(mfrow = c(3,1))
hist(deseq.alt[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : alt test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(deseq.null[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : null test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(pval[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : pvalue ", length(pval)-length(del.ix.deseq)), xlim=c(xmin, xmax), xlab ="p=value")




dev.off()






#################################################
## QQ plot of p-values (statistics) from null and alternative in DESeq
#################################################

setwd("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/fig/DESeqdebug/")
pdf("QQ.statistic.DESeq.debug.pdf")

#### filter.cut = 0

filter.cut = 0

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.txt"))[,1])
 
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))



numTests = length(deseq.alt[-del.ix.deseq])
X <- -log(((1:numTests)/numTests), 10)

num = 50000
xval = X[1:num]
yval.alt = sort(-log(deseq.alt[-del.ix.deseq],10), decreasing=T)[1:num]
yval.null = sort(-log(deseq.null[-del.ix.deseq],10), decreasing=T)[1:num]

xymax = max(xval, yval.null, yval.alt)

plot(xval, yval.alt, xlim=c(0,xymax), ylim=c(0, xymax), cex = 0.3, ylab = "log10 quantiles of p-value distribution",xlab = "log10 quantiles of uniform distribution", main=paste0("filter ", filter.cut), col="red")
abline(0,1)
points(xval, yval.null, col="blue", cex = 0.3)
legend(0,xymax, c("alt",  "null"), col = c("red", "blue"),pch = 19,cex = 1, text.col = "black",merge = FALSE, bg = "white")




#### filter.cut = 10

filter.cut = 10

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.", filter.cut, ".txt"))[,1])
  
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))



numTests = length(deseq.alt[-del.ix.deseq])
X <- -log(((1:numTests)/numTests), 10)

num = 50000
xval = X[1:num]
yval.alt = sort(-log(deseq.alt[-del.ix.deseq],10), decreasing=T)[1:num]
yval.null = sort(-log(deseq.null[-del.ix.deseq],10), decreasing=T)[1:num]

xymax = max(xval, yval.null, yval.alt)

plot(xval, yval.alt, xlim=c(0,xymax), ylim=c(0, xymax), cex = 0.3, ylab = "log10 quantiles of p-value distribution",xlab = "log10 quantiles of uniform distribution", main=paste0("filter ", filter.cut), col="red")
abline(0,1)
points(xval, yval.null, col="blue", cex = 0.3)
legend(0,xymax, c("alt",  "null"), col = c("red", "blue"),pch = 19,cex = 1, text.col = "black",merge = FALSE, bg = "white")






#### filter.cut = 20

filter.cut = 20

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.", filter.cut, ".txt"))[,1])
  
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))



numTests = length(deseq.alt[-del.ix.deseq])
X <- -log(((1:numTests)/numTests), 10)

num = 50000
xval = X[1:num]
yval.alt = sort(-log(deseq.alt[-del.ix.deseq],10), decreasing=T)[1:num]
yval.null = sort(-log(deseq.null[-del.ix.deseq],10), decreasing=T)[1:num]

xymax = max(xval, yval.null, yval.alt)

plot(xval, yval.alt, xlim=c(0,xymax), ylim=c(0, xymax), cex = 0.3, ylab = "log10 quantiles of p-value distribution",xlab = "log10 quantiles of uniform distribution", main=paste0("filter ", filter.cut), col="red")
abline(0,1)
points(xval, yval.null, col="blue", cex = 0.3)
legend(0,xymax, c("alt",  "null"), col = c("red", "blue"),pch = 19,cex = 1, text.col = "black",merge = FALSE, bg = "white")






#### filter.cut = 30

filter.cut = 30

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.", filter.cut, ".txt"))[,1])
  
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))



numTests = length(deseq.alt[-del.ix.deseq])
X <- -log(((1:numTests)/numTests), 10)

num = 50000
xval = X[1:num]
yval.alt = sort(-log(deseq.alt[-del.ix.deseq],10), decreasing=T)[1:num]
yval.null = sort(-log(deseq.null[-del.ix.deseq],10), decreasing=T)[1:num]

xymax = max(xval, yval.null, yval.alt)

plot(xval, yval.alt, xlim=c(0,xymax), ylim=c(0, xymax), cex = 0.3, ylab = "log10 quantiles of p-value distribution",xlab = "log10 quantiles of uniform distribution", main=paste0("filter ", filter.cut), col="red")
abline(0,1)
points(xval, yval.null, col="blue", cex = 0.3)
legend(0,xymax, c("alt",  "null"), col = c("red", "blue"),pch = 19,cex = 1, text.col = "black",merge = FALSE, bg = "white")






#### filter.cut = 60

filter.cut = 60

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.", filter.cut, ".txt"))[,1])
  
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))



numTests = length(deseq.alt[-del.ix.deseq])
X <- -log(((1:numTests)/numTests), 10)

num = 50000
xval = X[1:num]
yval.alt = sort(-log(deseq.alt[-del.ix.deseq],10), decreasing=T)[1:num]
yval.null = sort(-log(deseq.null[-del.ix.deseq],10), decreasing=T)[1:num]

xymax = max(xval, yval.null, yval.alt)

plot(xval, yval.alt, xlim=c(0,xymax), ylim=c(0, xymax), cex = 0.3, ylab = "log10 quantiles of p-value distribution",xlab = "log10 quantiles of uniform distribution", main=paste0("filter ", filter.cut), col="red")
abline(0,1)
points(xval, yval.null, col="blue", cex = 0.3)
legend(0,xymax, c("alt",  "null"), col = c("red", "blue"),pch = 19,cex = 1, text.col = "black",merge = FALSE, bg = "white")








dev.off()







#####################
### Storey FDR 
#####################

#install.packages("qvalue")
library("qvalue")

input.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/pval.ms.wave.", all.name, ".Robj")
load(input.path)
## pval.ms pval.wave


setwd("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/fig/DESeqdebug/")
pdf("discovery.FDR.both.pdf")

## filter.cut : 0
filter.cut = 0
pval.deseq = pval.deseq.3

del.ix.deseq = union(union(which(is.na(pval.deseq) == TRUE), which(is.na(pval.ms)==TRUE)), which(is.na(pval.wave)==TRUE))

qval.wave = qvalue(pval.wave[-del.ix.deseq])
qval.ms = qvalue(pval.ms[-del.ix.deseq])
qval.deseq = qvalue(pval.deseq[-del.ix.deseq])


qval.wave$pi0
qval.ms$pi0
qval.deseq$pi0

## 0.7818706
## 0.883166
## 0.9528259


alpha.list = seq(0.01, 0.2, by=0.01)
length(alpha.list)
# 20
num.wave = num.ms = num.deseq = rep(NA, length(alpha.list))
for(i in 1:length(alpha.list)){
    num.wave[i] = sum(qval.wave$qvalues < alpha.list[i])
    num.ms[i] = sum(qval.ms$qvalues < alpha.list[i])
    num.deseq[i] = sum(qval.deseq$qvalues < alpha.list[i])
}


ymax = max(num.wave, num.ms, num.deseq)
ymin = 0
plot(alpha.list, num.ms, ylim=c(ymin,ymax), col="red", type = "l", main=paste0("filter : ", filter.cut), xlab = "FDR", ylab="number of significant tests")
points(alpha.list, num.wave, ylim=c(ymin,ymax), col="blue", type="l")
points(alpha.list, num.deseq, ylim=c(ymin,ymax), col="darkgreen", type="l")
legend(0.01, ymax*0.99, c("multiseq", "WaveQTL", "DESeq") , col = c("red", "blue", "darkgreen"), lty = c(1,1,1), text.col = "black",merge = FALSE, bg = "white")



## filter.cut : 10
filter.cut = 10
pval.deseq = pval.deseq.3.10

del.ix.deseq = union(union(which(is.na(pval.deseq) == TRUE), which(is.na(pval.ms)==TRUE)), which(is.na(pval.wave)==TRUE))

qval.wave = qvalue(pval.wave[-del.ix.deseq])
qval.ms = qvalue(pval.ms[-del.ix.deseq])
qval.deseq = qvalue(pval.deseq[-del.ix.deseq])


qval.wave$pi0
qval.ms$pi0
qval.deseq$pi0

## 0.7398409
## 0.8730911
## 0.8939121

alpha.list = seq(0.01, 0.2, by=0.01)
length(alpha.list)
# 20
num.wave = num.ms = num.deseq = rep(NA, length(alpha.list))
for(i in 1:length(alpha.list)){
    num.wave[i] = sum(qval.wave$qvalues < alpha.list[i])
    num.ms[i] = sum(qval.ms$qvalues < alpha.list[i])
    num.deseq[i] = sum(qval.deseq$qvalues < alpha.list[i])
}


ymax = max(num.wave, num.ms, num.deseq)
ymin = 0
plot(alpha.list, num.ms, ylim=c(ymin,ymax), col="red", type = "l", main=paste0("filter : ", filter.cut), xlab = "FDR", ylab="number of significant tests")
points(alpha.list, num.wave, ylim=c(ymin,ymax), col="blue", type="l")
points(alpha.list, num.deseq, ylim=c(ymin,ymax), col="darkgreen", type="l")
legend(0.01, ymax*0.99, c("multiseq", "WaveQTL", "DESeq") , col = c("red", "blue", "darkgreen"), lty = c(1,1,1), text.col = "black",merge = FALSE, bg = "white")




## filter.cut : 20
filter.cut = 20
pval.deseq = pval.deseq.3.20

del.ix.deseq = union(union(which(is.na(pval.deseq) == TRUE), which(is.na(pval.ms)==TRUE)), which(is.na(pval.wave)==TRUE))

qval.wave = qvalue(pval.wave[-del.ix.deseq])
qval.ms = qvalue(pval.ms[-del.ix.deseq])
qval.deseq = qvalue(pval.deseq[-del.ix.deseq])


qval.wave$pi0
qval.ms$pi0
qval.deseq$pi0

## 0.7061634
## 0.8600601
## 0.7604863


alpha.list = seq(0.01, 0.2, by=0.01)
length(alpha.list)
# 20
num.wave = num.ms = num.deseq = rep(NA, length(alpha.list))
for(i in 1:length(alpha.list)){
    num.wave[i] = sum(qval.wave$qvalues < alpha.list[i])
    num.ms[i] = sum(qval.ms$qvalues < alpha.list[i])
    num.deseq[i] = sum(qval.deseq$qvalues < alpha.list[i])
}


ymax = max(num.wave, num.ms, num.deseq)
ymin = 0
plot(alpha.list, num.ms, ylim=c(ymin,ymax), col="red", type = "l", main=paste0("filter : ", filter.cut), xlab = "FDR", ylab="number of significant tests")
points(alpha.list, num.wave, ylim=c(ymin,ymax), col="blue", type="l")
points(alpha.list, num.deseq, ylim=c(ymin,ymax), col="darkgreen", type="l")
legend(0.01, ymax*0.99, c("multiseq", "WaveQTL", "DESeq") , col = c("red", "blue", "darkgreen"), lty = c(1,1,1), text.col = "black",merge = FALSE, bg = "white")






## filter.cut : 30
filter.cut = 30
pval.deseq = pval.deseq.3.30

del.ix.deseq = union(union(which(is.na(pval.deseq) == TRUE), which(is.na(pval.ms)==TRUE)), which(is.na(pval.wave)==TRUE))

qval.wave = qvalue(pval.wave[-del.ix.deseq])
qval.ms = qvalue(pval.ms[-del.ix.deseq])
qval.deseq = qvalue(pval.deseq[-del.ix.deseq])


qval.wave$pi0
qval.ms$pi0
qval.deseq$pi0

## 0.6372468
## 0.8163906
## 0.650188


alpha.list = seq(0.01, 0.2, by=0.01)
length(alpha.list)
# 20
num.wave = num.ms = num.deseq = rep(NA, length(alpha.list))
for(i in 1:length(alpha.list)){
    num.wave[i] = sum(qval.wave$qvalues < alpha.list[i])
    num.ms[i] = sum(qval.ms$qvalues < alpha.list[i])
    num.deseq[i] = sum(qval.deseq$qvalues < alpha.list[i])
}


ymax = max(num.wave, num.ms, num.deseq)
ymin = 0
plot(alpha.list, num.ms, ylim=c(ymin,ymax), col="red", type = "l", main=paste0("filter : ", filter.cut), xlab = "FDR", ylab="number of significant tests")
points(alpha.list, num.wave, ylim=c(ymin,ymax), col="blue", type="l")
points(alpha.list, num.deseq, ylim=c(ymin,ymax), col="darkgreen", type="l")
legend(0.01, ymax*0.99, c("multiseq", "WaveQTL", "DESeq") , col = c("red", "blue", "darkgreen"), lty = c(1,1,1), text.col = "black",merge = FALSE, bg = "white")






## filter.cut : 60
filter.cut = 60
pval.deseq = pval.deseq.3.60

del.ix.deseq = union(union(which(is.na(pval.deseq) == TRUE), which(is.na(pval.ms)==TRUE)), which(is.na(pval.wave)==TRUE))

qval.wave = qvalue(pval.wave[-del.ix.deseq])
qval.ms = qvalue(pval.ms[-del.ix.deseq])
qval.deseq = qvalue(pval.deseq[-del.ix.deseq])


qval.wave$pi0
qval.ms$pi0
qval.deseq$pi0

## 0.362377
## 0.5406936
## 0.7918262


alpha.list = seq(0.01, 0.2, by=0.01)
length(alpha.list)
# 20
num.wave = num.ms = num.deseq = rep(NA, length(alpha.list))
for(i in 1:length(alpha.list)){
    num.wave[i] = sum(qval.wave$qvalues < alpha.list[i])
    num.ms[i] = sum(qval.ms$qvalues < alpha.list[i])
    num.deseq[i] = sum(qval.deseq$qvalues < alpha.list[i])
}


ymax = max(num.wave, num.ms, num.deseq)
ymin = 0
plot(alpha.list, num.ms, ylim=c(ymin,ymax), col="red", type = "l", main=paste0("filter : ", filter.cut), xlab = "FDR", ylab="number of significant tests")
points(alpha.list, num.wave, ylim=c(ymin,ymax), col="blue", type="l")
points(alpha.list, num.deseq, ylim=c(ymin,ymax), col="darkgreen", type="l")
legend(0.01, ymax*0.99, c("multiseq", "WaveQTL", "DESeq") , col = c("red", "blue", "darkgreen"), lty = c(1,1,1), text.col = "black",merge = FALSE, bg = "white")


dev.off()




























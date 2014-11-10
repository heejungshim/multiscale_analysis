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


out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.noC.pval.discrete.Robj")

load(out.path)
## pval.deseq.3.0
## pval.deseq.3.60
## pval.deseq.3.30
## pval.deseq.3.20
## pval.deseq.3.10





#################################################
## make a histogram of statistics from null and alternative and p-values for DESeq
#################################################

setwd("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/fig/DESeqdebug/")
pdf("hist.statistic.pval.DESeq2.noC.debug.discrete.pdf")

#### filter.cut = 0

filter.cut = 0
pval = pval.deseq.3.0

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.noC.min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/DESeq2.noC.min.pval.", filter.cut, ".txt"))[,1])

 
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

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.noC.min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/DESeq2.noC.min.pval.", filter.cut, ".txt"))[,1])
 
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

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.noC.min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/DESeq2.noC.min.pval.", filter.cut, ".txt"))[,1])
 
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

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.noC.min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/DESeq2.noC.min.pval.", filter.cut, ".txt"))[,1])
 
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

deseq.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.noC.min.pval.", filter.cut, ".txt"))[,1])
deseq.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/DESeq2.noC.min.pval.", filter.cut, ".txt"))[,1])
 
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))

xmax = 1
xmin = 0

par(mfrow = c(3,1))
hist(deseq.alt[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : alt test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(deseq.null[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : null test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(pval[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : pvalue ", length(pval)-length(del.ix.deseq)), xlim=c(xmin, xmax), xlab ="p=value")




dev.off()








#################################################
## make a histogram of statistics from null and alternative and p-values for DESeq (pooling)
#################################################

setwd("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/fig/DESeqdebug/")
pdf("hist.statistic.pval.DESeq2.noC.debug.discrete.pooling.pdf")



#### filter.cut = 0

filter.cut = 0


deseq.alt = as.numeric(as.matrix(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.noC.all.pval.", filter.cut, ".txt"))))
deseq.null = as.numeric(as.matrix(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/DESeq2.noC.all.pval.", filter.cut, ".txt"))))

del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))


xmax = 1
xmin = 0

par(mfrow = c(2,1))
hist(deseq.alt[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : alt test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(deseq.null[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : null test statistic"), xlim=c(xmin, xmax), xlab ="p-value")

#### filter.cut = 10

filter.cut = 10

deseq.alt = as.numeric(as.matrix(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.noC.all.pval.", filter.cut, ".txt"))))
deseq.null = as.numeric(as.matrix(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/DESeq2.noC.all.pval.", filter.cut, ".txt"))))
  
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))


xmax = 1
xmin = 0

par(mfrow = c(2,1))
hist(deseq.alt[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : alt test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(deseq.null[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : null test statistic"), xlim=c(xmin, xmax), xlab ="p-value")


#### filter.cut = 20

filter.cut = 20

deseq.alt = as.numeric(as.matrix(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.noC.all.pval.", filter.cut, ".txt"))))
deseq.null = as.numeric(as.matrix(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/DESeq2.noC.all.pval.", filter.cut, ".txt"))))
  
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))


xmax = 1
xmin = 0

par(mfrow = c(2,1))
hist(deseq.alt[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : alt test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(deseq.null[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : null test statistic"), xlim=c(xmin, xmax), xlab ="p-value")

#### filter.cut = 30

filter.cut = 30

deseq.alt = as.numeric(as.matrix(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.noC.all.pval.", filter.cut, ".txt"))))
deseq.null = as.numeric(as.matrix(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/DESeq2.noC.all.pval.", filter.cut, ".txt"))))
  
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))




xmax = 1
xmin = 0

par(mfrow = c(2,1))
hist(deseq.alt[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : alt test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(deseq.null[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : null test statistic"), xlim=c(xmin, xmax), xlab ="p-value")





#### filter.cut = 60

filter.cut = 60

deseq.alt = as.numeric(as.matrix(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.noC.all.pval.", filter.cut, ".txt"))))
deseq.null = as.numeric(as.matrix(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/DESeq2.noC.all.pval.", filter.cut, ".txt"))))
  
del.ix.deseq = union(which(is.na(deseq.alt)==TRUE), which(is.na(deseq.null)==TRUE))



xmax = 1
xmin = 0

par(mfrow = c(2,1))
hist(deseq.alt[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : alt test statistic"), xlim=c(xmin, xmax), xlab ="p-value")
hist(deseq.null[-del.ix.deseq], breaks = 200, main=paste0(filter.cut, " : null test statistic"), xlim=c(xmin, xmax), xlab ="p-value")



dev.off()







#####################
### Storey FDR 
#####################

#install.packages("qvalue")
library("qvalue")


siteSize=2048
treatment='Copper'
strand='both'
window.size=300
numSam = 6
all.name = paste0(treatment,".", siteSize, ".", strand)



out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.noC.pval.discrete.Robj")

load(out.path)


input.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/pval.ms.wave.", all.name, ".Robj")
load(input.path)
## pval.ms pval.wave


setwd("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/fig/DESeqdebug/")
pdf("FDR.DESeq2.noC.discrete.pdf")

## filter.cut : 0
filter.cut = 0
pval.deseq = pval.deseq.3.0

del.ix.deseq = union(union(which(is.na(pval.deseq) == TRUE), which(is.na(pval.ms)==TRUE)), which(is.na(pval.wave)==TRUE))

qval.wave = qvalue(pval.wave[-del.ix.deseq])
qval.ms = qvalue(pval.ms[-del.ix.deseq])
qval.deseq = qvalue(pval.deseq[-del.ix.deseq])


qval.wave$pi0
qval.ms$pi0
qval.deseq$pi0

## 0.7818706
## 0.883166
## 0.9422277


alpha.list = seq(0, 0.2, by=0.01)
length(alpha.list)
# 20
num.wave = num.ms = num.deseq = rep(NA, length(alpha.list))
for(i in 1:length(alpha.list)){
    num.wave[i] = sum(qval.wave$qvalues < alpha.list[i])
    num.ms[i] = sum(qval.ms$qvalues < alpha.list[i])
    num.deseq[i] = sum(qval.deseq$qvalues < alpha.list[i])
}




num.wave = num.wave/length(qval.wave$qvalues)
num.ms = num.ms/length(qval.wave$qvalues)
num.deseq = num.deseq/length(qval.wave$qvalues)




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
## 0.884924

alpha.list = seq(0, 0.2, by=0.01)
length(alpha.list)
# 20
num.wave = num.ms = num.deseq = rep(NA, length(alpha.list))
for(i in 1:length(alpha.list)){
    num.wave[i] = sum(qval.wave$qvalues < alpha.list[i])
    num.ms[i] = sum(qval.ms$qvalues < alpha.list[i])
    num.deseq[i] = sum(qval.deseq$qvalues < alpha.list[i])
}



num.wave = num.wave/length(qval.wave$qvalues)
num.ms = num.ms/length(qval.wave$qvalues)
num.deseq = num.deseq/length(qval.wave$qvalues)




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
## 0.7507062


alpha.list = seq(0, 0.2, by=0.01)
length(alpha.list)
# 20
num.wave = num.ms = num.deseq = rep(NA, length(alpha.list))
for(i in 1:length(alpha.list)){
    num.wave[i] = sum(qval.wave$qvalues < alpha.list[i])
    num.ms[i] = sum(qval.ms$qvalues < alpha.list[i])
    num.deseq[i] = sum(qval.deseq$qvalues < alpha.list[i])
}




num.wave = num.wave/length(qval.wave$qvalues)
num.ms = num.ms/length(qval.wave$qvalues)
num.deseq = num.deseq/length(qval.wave$qvalues)



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
## 0.6760918


alpha.list = seq(0, 0.2, by=0.01)
length(alpha.list)
# 20
num.wave = num.ms = num.deseq = rep(NA, length(alpha.list))
for(i in 1:length(alpha.list)){
    num.wave[i] = sum(qval.wave$qvalues < alpha.list[i])
    num.ms[i] = sum(qval.ms$qvalues < alpha.list[i])
    num.deseq[i] = sum(qval.deseq$qvalues < alpha.list[i])
}




num.wave = num.wave/length(qval.wave$qvalues)
num.ms = num.ms/length(qval.wave$qvalues)
num.deseq = num.deseq/length(qval.wave$qvalues)



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
## 0.8301283


alpha.list = seq(0, 0.2, by=0.01)
length(alpha.list)
# 20
num.wave = num.ms = num.deseq = rep(NA, length(alpha.list))
for(i in 1:length(alpha.list)){
    num.wave[i] = sum(qval.wave$qvalues < alpha.list[i])
    num.ms[i] = sum(qval.ms$qvalues < alpha.list[i])
    num.deseq[i] = sum(qval.deseq$qvalues < alpha.list[i])
}




num.wave = num.wave/length(qval.wave$qvalues)
num.ms = num.ms/length(qval.wave$qvalues)
num.deseq = num.deseq/length(qval.wave$qvalues)



ymax = max(num.wave, num.ms, num.deseq)
ymin = 0
plot(alpha.list, num.ms, ylim=c(ymin,ymax), col="red", type = "l", main=paste0("filter : ", filter.cut), xlab = "FDR", ylab="number of significant tests")
points(alpha.list, num.wave, ylim=c(ymin,ymax), col="blue", type="l")
points(alpha.list, num.deseq, ylim=c(ymin,ymax), col="darkgreen", type="l")
legend(0.01, ymax*0.99, c("multiseq", "WaveQTL", "DESeq") , col = c("red", "blue", "darkgreen"), lty = c(1,1,1), text.col = "black",merge = FALSE, bg = "white")


dev.off()





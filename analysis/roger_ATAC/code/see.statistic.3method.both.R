#!/usr/bin/env Rscript

## Aim : This file contains Rscript to look at statistic patterns from null and alternative (three methods)
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+

siteSize=2048
treatment='Copper'
strand='both'

## directory name 
alt.name=paste0(treatment,".", siteSize, ".", strand,  ".alt")
null.name=paste0(treatment,".", siteSize, ".", strand, ".null")
all.name = paste0(treatment,".", siteSize, ".", strand)


## read p-value from DESeq

deseq.100.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".alt.run/output/min.pval.txt"))[,1])

length(deseq.100.alt)

deseq.100.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".null.run/output/min.pval.txt"))[,1])

length(deseq.100.null)



deseq.300.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/min.pval.txt"))[,1])

length(deseq.300.alt)

deseq.300.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/min.pval.txt"))[,1])

length(deseq.300.null)





## read logLR from wavelet
## alt
load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/wave/", alt.name, ".run/sum/", alt.name, ".Robj"))
wave.alt = unlist(logLR_list)
wave.done.alt = unlist(done_list)

## null
load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/wave/", null.name, ".run/sum/", null.name, ".Robj"))
wave.null = unlist(logLR_list)
wave.done.null = unlist(done_list)



## read logLR from multiseq
## alt
load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/multiscale/", alt.name, ".sum/", alt.name, ".Robj"))
ms.alt = unlist(logLR_list)
ms.done.alt = unlist(done_list)

## null
load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/multiscale/", null.name, ".sum/", null.name, ".Robj"))
ms.null = unlist(logLR_list)
ms.done.null = unlist(done_list)




#########################
# see invalid results
#########################

length(wave.alt)
length(wave.null)
length(wave.done.alt)
length(wave.done.null)
min(wave.alt)
max(wave.alt)
min(wave.null)
max(wave.null)
sum(is.na(wave.done.alt))
sum(wave.done.alt)
sum(is.na(wave.done.null))
sum(wave.done.null)


#> length(wave.alt)
#[1] 197969
#> length(wave.null)
#[1] 197969
#> length(wave.done.alt)
#[1] 197969
#> length(wave.done.null)
#[1] 197969
#> min(wave.alt)
#[1] 0
#> max(wave.alt)
#[1] 52.6089
#> min(wave.null)
#[1] 0
#> max(wave.null)
#[1] 11.38442
#> sum(is.na(wave.done.alt))
#[1] 0
#> sum(wave.done.alt)
#[1] 197969
#> sum(is.na(wave.done.null))
#[1] 0
#> sum(wave.done.null)
#[1] 197969



length(ms.alt)
length(ms.null)
length(ms.done.alt)
length(ms.done.null)
min(ms.alt)
max(ms.alt)
min(ms.null)
max(ms.null)
sum(is.na(ms.done.alt))
sum(ms.done.alt)
sum(is.na(ms.done.null))
sum(ms.done.null)



#> length(ms.alt)
#[1] 197969
#> length(ms.null)
#[1] 197969
#> length(ms.done.alt)
#[1] 197969
#> length(ms.done.null)
#[1] 197969
#> min(ms.alt)
#[1] NA
#> max(ms.alt)
#[1] NA
#> min(ms.null)
#[1] NA
#> max(ms.null)
#[1] NA
#> sum(is.na(ms.done.alt))
#[1] 0
#> sum(ms.done.alt)
#[1] 197968
#> sum(is.na(ms.done.null))
#[1] 0
#> sum(ms.done.null)
#[1] 197965


min(ms.alt, na.rm = TRUE)
max(ms.alt, na.rm = TRUE)
min(ms.null, na.rm = TRUE)
max(ms.null, na.rm = TRUE)

#> min(ms.alt, na.rm = TRUE)
#[1] 0
#> max(ms.alt, na.rm = TRUE)
#[1] 514.096
#> min(ms.null, na.rm = TRUE)
#[1] 0
#> max(ms.null, na.rm = TRUE)
#[1] 18.85504

length(deseq.100.alt)
length(deseq.300.alt)
length(deseq.100.null)
length(deseq.300.null)
sum(is.infinite(deseq.100.alt))
sum(is.infinite(deseq.300.alt))
sum(is.infinite(deseq.100.null))
sum(is.infinite(deseq.300.null))


#> sum(is.infinite(deseq.100.alt))
#[1] 11447
#> sum(is.infinite(deseq.300.alt))
#[1] 11447
#> sum(is.infinite(deseq.100.null))
#[1] 11462
#> sum(is.infinite(deseq.300.null))
#[1] 11462



deseq.100.alt.inf = which(is.infinite(deseq.100.alt)==TRUE)
deseq.300.alt.inf = which(is.infinite(deseq.300.alt)==TRUE)
deseq.100.null.inf = which(is.infinite(deseq.100.null)==TRUE)
deseq.300.null.inf = which(is.infinite(deseq.300.null)==TRUE)

sum(deseq.100.alt.inf != deseq.300.alt.inf)
# 0
sum(deseq.100.null.inf != deseq.300.null.inf)
# 0

sum(wave.alt[deseq.100.alt.inf])
# 0 
sum(wave.null[deseq.100.null.inf])
# 0

sum(ms.alt[deseq.100.alt.inf])
# 0 
sum(ms.null[deseq.100.null.inf])
# 0






del.ix = union(union(union(deseq.100.alt.inf, deseq.100.null.inf), which(is.na(ms.alt)==TRUE)), which(is.na(ms.null)==TRUE))
which(is.na(ms.alt)==TRUE)
which(is.na(ms.null)==TRUE)

#> which(is.na(ms.alt)==TRUE)
#[1] 35550
#> which(is.na(ms.null)==TRUE)
#[1]  61034 126327 126328 133432





#######################################
# select sites with at least one read
#######################################
wave.a = wave.alt[-del.ix]
wave.n = wave.null[-del.ix]
ms.a = ms.alt[-del.ix]
ms.n = ms.null[-del.ix]
deseq.100.a = deseq.100.alt[-del.ix]
deseq.300.a = deseq.300.alt[-del.ix]
deseq.100.n = deseq.100.null[-del.ix]
deseq.300.n = deseq.300.null[-del.ix]

    
#########################
# Make a histogram
#########################
 
setwd("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/code/")
pdf("hist.3.methods.both.pdf", width = 7, height = 5)

## Wave
xmax = max(wave.a, wave.n)
xmin = min(wave.a, wave.n)
wave.a.hist = hist(wave.a, plot=FALSE, breaks = 200)
wave.n.hist = hist(wave.n, plot=FALSE, breaks = 200)
ymax = max(wave.a.hist$counts, wave.n.hist$counts)
ymin = min(wave.a.hist$counts, wave.n.hist$counts)

par(mfrow = c(2,1))
hist(wave.n, breaks = 200, main="wave null", xlim=c(xmin, xmax), ylim=c(ymin, ymax))
hist(wave.a, breaks = 200, main="wave alt", xlim=c(xmin, xmax), ylim=c(ymin, ymax))


## ms
xmax = max(ms.a, ms.n)
xmin = min(ms.a, ms.n)
ms.a.hist = hist(ms.a, plot=FALSE, breaks = 200)
ms.n.hist = hist(ms.n, plot=FALSE, breaks = 200)
ymax = max(ms.a.hist$counts, ms.n.hist$counts)
ymin = min(ms.a.hist$counts, ms.n.hist$counts)

par(mfrow = c(2,1))
hist(ms.n, breaks = 200, main="multiseq null", xlim=c(xmin, xmax), ylim=c(ymin, ymax))
hist(ms.a, breaks = 200, main="multiseq alt", xlim=c(xmin, xmax), ylim=c(ymin, ymax))



## DEseq 100
xmax = 1
xmin = 0
deseq.100.a.hist = hist(deseq.100.a, plot=FALSE, breaks = 200)
deseq.100.n.hist = hist(deseq.100.n, plot=FALSE, breaks = 200)
ymax = max(deseq.100.a.hist$counts, deseq.100.n.hist$counts)
ymin = min(deseq.100.a.hist$counts, deseq.100.n.hist$counts)


par(mfrow = c(2,1))
hist(deseq.100.n, breaks = 200, main="DESeq 100 null", xlim=c(xmin, xmax), ylim=c(ymin, ymax))
hist(deseq.100.a, breaks = 200, main="DESeq 100 alt", xlim=c(xmin, xmax), ylim=c(ymin, ymax))



## DEseq 300
xmax = 1
xmin = 0
deseq.300.a.hist = hist(deseq.300.a, plot=FALSE, breaks = 200)
deseq.300.n.hist = hist(deseq.300.n, plot=FALSE, breaks = 200)
ymax = max(deseq.300.a.hist$counts, deseq.300.n.hist$counts)
ymin = min(deseq.300.a.hist$counts, deseq.300.n.hist$counts)


par(mfrow = c(2,1))
hist(deseq.300.n, breaks = 200, main="DESeq 300 null", xlim=c(xmin, xmax), ylim=c(ymin, ymax))
hist(deseq.300.a, breaks = 200, main="DESeq 300 alt", xlim=c(xmin, xmax), ylim=c(ymin, ymax))


dev.off()



#################################
# put randomly selected (from unif[-0.001, 0] logLR for 0
#################################
set.seed(1)
ms.new.a = ms.a
wh = which(ms.a == 0)
length(wh)
ms.new.a[wh] = runif(length(wh), -0.001,0)

ms.new.n = ms.n
wh = which(ms.n == 0)
length(wh)
ms.new.n[wh] = runif(length(wh), -0.001,0)

wave.new.a = wave.a
wh = which(wave.a == 0)
length(wh)
wave.new.a[wh] = runif(length(wh), -0.001,0)

wave.new.n = wave.n
wh = which(wave.n == 0)
length(wh)
wave.new.n[wh] = runif(length(wh), -0.001,0)




##################################################################
# get p-value from empirical null distribution of test statistic
##################################################################


get.pval.from.empirical.null.dist <- function(statistic.null, statistic.alt, big.sig = TRUE){
    
    numNulltests = length(statistic.null)
    if(big.sig){
        numSig = sapply(statistic.alt, function(x, statistic.null){ return(sum(statistic.null >= x))}, statistic.null = statistic.null)
    }else{
        numSig = sapply(statistic.alt, function(x, statistic.null){ return(sum(statistic.null <= x))}, statistic.null = statistic.null)
    }
    pval.list = sapply(numSig, function(x, numNulltests){ return(runif(1, (x+1)/(numNulltests+2), (x+1)/(numNulltests+1)))}, numNulltests = numNulltests)

    return(pval.list)
}

pval.wave = get.pval.from.empirical.null.dist(statistic.null = wave.n, statistic.alt = wave.a)

pval.ms = get.pval.from.empirical.null.dist(statistic.null = ms.n, statistic.alt = ms.a)

pval.wave.new = get.pval.from.empirical.null.dist(statistic.null = wave.new.n, statistic.alt = wave.new.a)

pval.ms.new = get.pval.from.empirical.null.dist(statistic.null = ms.new.n, statistic.alt = ms.new.a)


pval.deseq.100 = get.pval.from.empirical.null.dist(statistic.null = deseq.100.n, statistic.alt = deseq.100.a, big.sig = FALSE)

pval.deseq.300 = get.pval.from.empirical.null.dist(statistic.null = deseq.300.n, statistic.alt = deseq.300.a, big.sig = FALSE)


save("pval.deseq.100", "pval.deseq.300", file =paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/deseq.", all.name, ".Robj"))

     
save("pval.wave", file =paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/wave.", all.name, ".Robj"))

save("pval.ms", file =paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/ms.", all.name, ".Robj"))

save("pval.wave.new", file =paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/wave.new.", all.name, ".Robj"))

save("pval.ms.new", file =paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/ms.new.", all.name, ".Robj"))



#################################################
## make a histogram of p-values, fdr, QQ-plot? 
#################################################

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/deseq.", all.name, ".Robj"))

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/wave.", all.name, ".Robj"))

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/ms.", all.name, ".Robj"))


setwd("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/code/")
pdf("hist.pval.3.methods.both.pdf")

xmax = 1
xmin = 0

par(mfrow = c(4,1))
hist(pval.wave, breaks = 200, main="wave", xlim=c(xmin, xmax))
hist(pval.ms, breaks = 200, main="multiseq", xlim=c(xmin, xmax))
hist(pval.deseq.100, breaks = 200, main="DESeq 100", xlim=c(xmin, xmax))
hist(pval.deseq.300, breaks = 200, main="DESeq 300", xlim=c(xmin, xmax))

dev.off()



#################################################
## make a histogram of p-values with new
#################################################

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/deseq.", all.name, ".Robj"))

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/wave.new.", all.name, ".Robj"))

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/ms.new.", all.name, ".Robj"))


setwd("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/code/")
pdf("hist.pval.3.methods.new.both.pdf")

xmax = 1
xmin = 0

par(mfrow = c(4,1))
hist(pval.wave.new, breaks = 200, main="wave", xlim=c(xmin, xmax))
hist(pval.ms.new, breaks = 200, main="multiseq", xlim=c(xmin, xmax))
hist(pval.deseq.100, breaks = 200, main="DESeq 100", xlim=c(xmin, xmax))
hist(pval.deseq.300, breaks = 200, main="DESeq 300", xlim=c(xmin, xmax))

dev.off()



#################################
# BH fdr
#################################

pval.bh.wave = p.adjust(pval.wave.new, method="BH", n= length(pval.wave.new))
pval.bh.ms = p.adjust(pval.ms.new, method="BH", n= length(pval.ms.new))
pval.bh.deseq.100 = p.adjust(pval.deseq.100, method="BH", n= length(pval.deseq.100))
pval.bh.deseq.300 = p.adjust(pval.deseq.300, method="BH", n= length(pval.deseq.300))



##################################################
# number of discovery for BH fdr
##################################################


alpha.list = seq(0.01, 0.2, by=0.01)
length(alpha.list)
# 20
num.wave = num.ms = num.deseq.100 = num.deseq.300 = rep(NA, length(alpha.list))
for(i in 1:length(alpha.list)){
    num.wave[i] = sum(pval.bh.wave < alpha.list[i])
    num.ms[i] = sum(pval.bh.ms < alpha.list[i])
    num.deseq.100[i] = sum(pval.bh.deseq.100 < alpha.list[i])
    num.deseq.300[i] = sum(pval.bh.deseq.300 < alpha.list[i])
}


setwd("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/code/")
pdf("discovery.BH.both.pdf")

ymax = max(num.wave, num.ms, num.deseq.100, num.deseq.300)
ymin = 0
plot(alpha.list, num.ms, ylim=c(ymin,ymax), col="red", type = "l")
points(alpha.list, num.wave, ylim=c(ymin,ymax), col="blue", type="l")
points(alpha.list, num.deseq.100, ylim=c(ymin,ymax), col="darkgreen", type="l")
points(alpha.list, num.deseq.300, ylim=c(ymin,ymax), col="darkgreen", type="l", lty="dashed")

dev.off()






###########################
# check intersection
###########################


alpha = 0.2

sum(pval.bh.wave< alpha)
sum(pval.bh.ms < alpha)
sum(pval.bh.deseq.300 < alpha)

sum((pval.bh.wave < alpha) & (pval.bh.ms < alpha))
sum((pval.bh.wave < alpha) & (pval.bh.deseq.300 < alpha))
sum((pval.bh.ms < alpha) & (pval.bh.deseq.300 < alpha))

sum((pval.bh.wave < alpha) & (pval.bh.ms < alpha))/sum(pval.bh.ms < alpha)
sum((pval.bh.wave < alpha) & (pval.bh.ms < alpha))/sum(pval.bh.wave < alpha)

sum((pval.bh.deseq.300 < alpha) & (pval.bh.ms < alpha))/sum(pval.bh.ms < alpha)
sum((pval.bh.deseq.300 < alpha) & (pval.bh.ms < alpha))/sum(pval.bh.deseq.300 < alpha)

sum((pval.bh.deseq.300 < alpha) & (pval.bh.wave < alpha))/sum(pval.bh.wave < alpha)
sum((pval.bh.deseq.300 < alpha) & (pval.bh.wave < alpha))/sum(pval.bh.deseq.300 < alpha)



#> sum(pval.bh.wave< alpha)
#[1] 275
#> sum(pval.bh.ms < alpha)
#[1] 2895
#> sum(pval.bh.deseq.300 < alpha)
#[1] 71
 
#> sum((pval.bh.wave < alpha) & (pval.bh.ms < alpha))
#[1] 204
#> sum((pval.bh.wave < alpha) & (pval.bh.deseq.300 < alpha))
#[1] 27
#> sum((pval.bh.ms < alpha) & (pval.bh.deseq.300 < alpha))
#[1] 60
 
#> sum((pval.bh.wave < alpha) & (pval.bh.ms < alpha))/sum(pval.bh.ms < alpha)
#[1] 0.07046632
#> sum((pval.bh.wave < alpha) & (pval.bh.ms < alpha))/sum(pval.bh.wave < alpha)
#[1] 0.7418182
 
#> sum((pval.bh.deseq.300 < alpha) & (pval.bh.ms < alpha))/sum(pval.bh.ms < alpha)
#[1] 0.02072539
#> sum((pval.bh.deseq.300 < alpha) & (pval.bh.ms < alpha))/sum(pval.bh.deseq.300 < alpha)
#[1] 0.8450704
 
#> sum((pval.bh.deseq.300 < alpha) & (pval.bh.wave < alpha))/sum(pval.bh.wave < alpha)
#[1] 0.09818182
#> sum((pval.bh.deseq.300 < alpha) & (pval.bh.wave < alpha))/sum(pval.bh.deseq.300 < alpha)
#[1] 0.3802817



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


alpha.list = seq(0.01, 0.2, by=0.01)
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
pdf("discovery.FDR.both.pdf")

ymax = max(num.wave, num.ms, num.deseq.100, num.deseq.300)
ymin = 0
plot(alpha.list, num.ms, ylim=c(ymin,ymax), col="red", type = "l")
points(alpha.list, num.wave, ylim=c(ymin,ymax), col="blue", type="l")
points(alpha.list, num.deseq.100, ylim=c(ymin,ymax), col="darkgreen", type="l")
points(alpha.list, num.deseq.300, ylim=c(ymin,ymax), col="darkgreen", type="l", lty="dashed")

dev.off()



###########################
# check intersection
###########################


alpha = 0.2

sum(qval.wave$qvalues < alpha)
sum(qval.ms$qvalues < alpha)
sum(qval.deseq.300$qvalues < alpha)

sum((qval.wave$qvalues < alpha) & (qval.ms$qvalues < alpha))
sum((qval.wave$qvalues < alpha) & (qval.deseq.300$qvalues < alpha))
sum((qval.ms$qvalues < alpha) & (qval.deseq.300$qvalues < alpha))

sum((qval.wave$qvalues < alpha) & (qval.ms$qvalues < alpha))/sum(qval.ms$qvalues < alpha)
sum((qval.wave$qvalues < alpha) & (qval.ms$qvalues < alpha))/sum(qval.wave$qvalues < alpha)

sum((qval.deseq.300$qvalues < alpha) & (qval.ms$qvalues < alpha))/sum(qval.ms$qvalues < alpha)
sum((qval.deseq.300$qvalues < alpha) & (qval.ms$qvalues < alpha))/sum(qval.deseq.300$qvalues < alpha)

sum((qval.deseq.300$qvalues < alpha) & (qval.wave$qvalues < alpha))/sum(qval.wave$qvalues < alpha)
sum((qval.deseq.300$qvalues < alpha) & (qval.wave$qvalues < alpha))/sum(qval.deseq.300$qvalues < alpha)


#> sum(qval.wave$qvalues < alpha)
#[1] 377
#> sum(qval.ms$qvalues < alpha)
#[1] 3477
#> sum(qval.deseq.300$qvalues < alpha)
#[1] 71
 
#> sum((qval.wave$qvalues < alpha) & (qval.ms$qvalues < alpha))
#[1] 301
#> sum((qval.wave$qvalues < alpha) & (qval.deseq.300$qvalues < alpha))
#[1] 30
#> sum((qval.ms$qvalues < alpha) & (qval.deseq.300$qvalues < alpha))
#[1] 61
 
#> sum((qval.wave$qvalues < alpha) & (qval.ms$qvalues < alpha))/sum(qval.ms$qvalues < alpha)
#[1] 0.08656888
#> sum((qval.wave$qvalues < alpha) & (qval.ms$qvalues < alpha))/sum(qval.wave$qvalues < alpha)
#[1] 0.7984085
 
#> sum((qval.deseq.300$qvalues < alpha) & (qval.ms$qvalues < alpha))/sum(qval.ms$qvalues < alpha)
#[1] 0.01754386
#> sum((qval.deseq.300$qvalues < alpha) & (qval.ms$qvalues < alpha))/sum(qval.deseq.300$qvalues < alpha)
#[1] 0.8591549
 
#> sum((qval.deseq.300$qvalues < alpha) & (qval.wave$qvalues < alpha))/sum(qval.wave$qvalues < alpha)
#[1] 0.0795756
#> sum((qval.deseq.300$qvalues < alpha) & (qval.wave$qvalues < alpha))/sum(qval.deseq.300$qvalues < alpha)
#[1] 0.4225352
 

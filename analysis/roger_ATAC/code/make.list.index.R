#!/usr/bin/env Rscript

## Aim : make a list of sites (index) to make effect size figure
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+


siteSize=2048
treatment='Copper'
strand='both'

## directory name
all.name = paste0(treatment,".", siteSize, ".", strand)

## read empirical p-values
## wave and ms
input.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/pval.ms.wave.", all.name, ".Robj")
load(input.path)
## pval.ms pval.wave


## DESeq2 300
out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.pval.discrete.Robj")
load(out.path)
## pval.deseq.3.0

## DESeq2 100 and 2048
input.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".alt.run/output/DESeq2.pval.discrete.Robj")
load(input.path)

## DESeq2 600
out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".alt.run/output/DESeq2.pval.discrete.",600,".Robj")

load(out.path)



#####################
## getting q-value
#####################

library("qvalue")


del.ix = NULL
del.ix = union(del.ix, which(is.na(pval.deseq.100.0)==TRUE))
del.ix = union(del.ix, which(is.na(pval.deseq.3.0)==TRUE))
del.ix = union(del.ix, which(is.na(pval.deseq.600.0)==TRUE))
del.ix = union(del.ix, which(is.na(pval.deseq.full.0)==TRUE))
del.ix = union(del.ix, which(is.na(pval.ms)==TRUE))
del.ix = union(del.ix, which(is.na(pval.wave)==TRUE))

length(del.ix)
## 11570

numTests = length(pval.ms)

qval.wave = qval.ms = qval.deseq.300 = qval.deseq.600 = qval.deseq.full = qval.deseq.100 = rep(NA, numTests)
index.f = (1:numTests)[-del.ix]

qval.wave[index.f] = qvalue(pval.wave[-del.ix])$qvalues
qval.ms[index.f]  = qvalue(pval.ms[-del.ix])$qvalues
qval.deseq.full[index.f] = qvalue(pval.deseq.full.0[-del.ix])$qvalues
qval.deseq.300[index.f] = qvalue(pval.deseq.3.0[-del.ix])$qvalues
qval.deseq.100[index.f] = qvalue(pval.deseq.100.0[-del.ix])$qvalues
qval.deseq.600[index.f] = qvalue(pval.deseq.600.0[-del.ix])$qvalues


output.path.all = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/tmp/", all.name) 

############################
## DESeq2 300 vs multiseq
############################
sum((qval.deseq.300 < 0.02) & (qval.ms < 0.01), na.rm = TRUE)
## 48
sum((qval.deseq.300 > 0.5) & (qval.ms < 0.01), na.rm = TRUE)
## 61
sum((qval.deseq.300 <  0.1) & (qval.ms > 0.4), na.rm = TRUE)
## 29

com.name = "DESeq2.3.ms"


ix = which((qval.deseq.300 < 0.02) & (qval.ms < 0.01))
path.each = paste0(output.path.all, ".", com.name, ".", "both.index")
cat(ix, file = path.each)

ix = which((qval.deseq.300 > 0.5) & (qval.ms < 0.01))
path.each = paste0(output.path.all, ".", com.name, ".", "ms.index")
cat(ix, file = path.each)

ix = which((qval.deseq.300 <  0.1) & (qval.ms > 0.4))
path.each = paste0(output.path.all, ".", com.name, ".", "DESeq2.3.index")
cat(ix, file = path.each)



############################
## DESeq2 600 vs multiseq
############################

sum((qval.deseq.600 < 0.05) & (qval.ms < 0.05), na.rm = TRUE)
## 67
sum((qval.deseq.600 > 0.6) & (qval.ms < 0.01), na.rm = TRUE)
## 47
sum((qval.deseq.600 <  0.1) & (qval.ms > 0.4), na.rm = TRUE)
## 11


com.name = "DESeq2.6.ms"

ix = which((qval.deseq.600 < 0.05) & (qval.ms < 0.05))
path.each = paste0(output.path.all, ".", com.name, ".", "both.index")
cat(ix, file = path.each)

ix = which((qval.deseq.600 > 0.6) & (qval.ms < 0.01))
path.each = paste0(output.path.all, ".", com.name, ".", "ms.index")
cat(ix, file = path.each)

ix = which((qval.deseq.600 <  0.1) & (qval.ms > 0.4))
path.each = paste0(output.path.all, ".", com.name, ".", "DESeq2.6.index")
cat(ix, file = path.each)




############################
## DESeq2 2048 vs multiseq
############################
sum((qval.deseq.full < 0.02) & (qval.ms < 0.01), na.rm = TRUE)
## 55
sum((qval.deseq.full > 0.5) & (qval.ms < 0.01), na.rm = TRUE)
## 29
sum((qval.deseq.full <  0.1) & (qval.ms > 0.4), na.rm = TRUE)
## 13

com.name = "DESeq2.full.ms"

ix = which((qval.deseq.full < 0.02) & (qval.ms < 0.01))
path.each = paste0(output.path.all, ".", com.name, ".", "both.index")
cat(ix, file = path.each)

ix = which((qval.deseq.full > 0.5) & (qval.ms < 0.01))
path.each = paste0(output.path.all, ".", com.name, ".", "ms.index")
cat(ix, file = path.each)

ix = which((qval.deseq.full <  0.1) & (qval.ms > 0.4))
path.each = paste0(output.path.all, ".", com.name, ".", "DESeq2.full.index")
cat(ix, file = path.each)



############################
## DESeq2 2048 vs DESeq2 600
############################
sum((qval.deseq.full < 0.02) & (qval.deseq.600 < 0.02), na.rm = TRUE)
## 24
sum((qval.deseq.full > 0.5) & (qval.deseq.600 < 0.05), na.rm = TRUE)
## 5
sum((qval.deseq.full <  0.05) & (qval.deseq.600 > 0.5), na.rm = TRUE)
## 6


com.name = "DESeq2.full.DESeq2.6"

ix = which((qval.deseq.full < 0.02) & (qval.deseq.600 < 0.02))
path.each = paste0(output.path.all, ".", com.name, ".", "both.index")
cat(ix, file = path.each)

ix = which((qval.deseq.full > 0.5) & (qval.deseq.600 < 0.05))
path.each = paste0(output.path.all, ".", com.name, ".", "DESeq2.6.index")
cat(ix, file = path.each)

ix = which((qval.deseq.full <  0.05) & (qval.deseq.600 > 0.5))
path.each = paste0(output.path.all, ".", com.name, ".", "DESeq2.full.index")
cat(ix, file = path.each)





############################
## DESeq2 2048 vs DESeq2 300
############################
sum((qval.deseq.full < 0.015) & (qval.deseq.300 < 0.015), na.rm = TRUE)
## 5
sum((qval.deseq.full > 0.5) & (qval.deseq.300 < 0.05), na.rm = TRUE)
## 34
sum((qval.deseq.full <  0.05) & (qval.deseq.300 > 0.5), na.rm = TRUE)
## 15



com.name = "DESeq2.full.DESeq2.3"

ix = which((qval.deseq.full < 0.015) & (qval.deseq.300 < 0.015))
path.each = paste0(output.path.all, ".", com.name, ".", "both.index")
cat(ix, file = path.each)

ix = which((qval.deseq.full > 0.5) & (qval.deseq.300 < 0.05))
path.each = paste0(output.path.all, ".", com.name, ".", "DESeq2.3.index")
cat(ix, file = path.each)


ix = which((qval.deseq.full <  0.05) & (qval.deseq.300 > 0.5))
path.each = paste0(output.path.all, ".", com.name, ".", "DESeq2.full.index")
cat(ix, file = path.each)




ix = which((qval.deseq.full < 0.02) & (qval.deseq.600 < 0.02))
path.each = paste0(output.path.all, ".", com.name, ".", "both.index")
cat(ix, file = path.each)





############################
## DESeq2 600 vs DESeq2 300
############################
sum((qval.deseq.600 < 0.02) & (qval.deseq.300 < 0.02), na.rm = TRUE)
## 29
sum((qval.deseq.600 > 0.5) & (qval.deseq.300 < 0.1), na.rm = TRUE)
## 24
sum((qval.deseq.600 <  0.1) & (qval.deseq.300 > 0.4), na.rm = TRUE)
## 2


com.name = "DESeq2.6.DESeq2.3"

ix = which((qval.deseq.600 < 0.02) & (qval.deseq.300 < 0.02))
path.each = paste0(output.path.all, ".", com.name, ".", "both.index")
cat(ix, file = path.each)

ix = which((qval.deseq.600 > 0.5) & (qval.deseq.300 < 0.1))
path.each = paste0(output.path.all, ".", com.name, ".", "DESeq2.3.index")
cat(ix, file = path.each)

ix = which((qval.deseq.600 <  0.1) & (qval.deseq.300 > 0.4))
path.each = paste0(output.path.all, ".", com.name, ".", "DESeq2.6.index")
cat(ix, file = path.each)


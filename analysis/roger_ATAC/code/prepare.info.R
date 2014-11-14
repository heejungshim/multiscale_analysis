#!/usr/bin/env Rscript

## Aim : This file contains Rscript to prepare information for figures
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+


## index is given!
## index.here :
## file.path : where (including file name) we want to save output!!!

##file.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/tmp/",  all.name, ".both.DESeq2.info")
##info.all = read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/tmp/",  all.name, ".both.DESeq.info"), header = TRUE)
##index.here = info.all$index


input.dir.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/tmp/"


index.files = c("Copper.2048.both.DESeq2.3.ms.both", "Copper.2048.both.DESeq2.3.ms.DESeq2.3", "Copper.2048.both.DESeq2.3.ms.ms", "Copper.2048.both.DESeq2.6.DESeq2.3.both", "Copper.2048.both.DESeq2.6.DESeq2.3.DESeq2.3", "Copper.2048.both.DESeq2.6.DESeq2.3.DESeq2.6", "Copper.2048.both.DESeq2.6.ms.both", "Copper.2048.both.DESeq2.6.ms.DESeq2.6", "Copper.2048.both.DESeq2.6.ms.ms", "Copper.2048.both.DESeq2.full.DESeq2.3.both", "Copper.2048.both.DESeq2.full.DESeq2.3.DESeq2.3", "Copper.2048.both.DESeq2.full.DESeq2.3.DESeq2.full", "Copper.2048.both.DESeq2.full.DESeq2.6.both", "Copper.2048.both.DESeq2.full.DESeq2.6.DESeq2.6", "Copper.2048.both.DESeq2.full.DESeq2.6.DESeq2.full", "Copper.2048.both.DESeq2.full.ms.both", "Copper.2048.both.DESeq2.full.ms.DESeq2.full", "Copper.2048.both.DESeq2.full.ms.ms")




for(sss in 1:length(index.files)){

  ##sss = 1
  index.here = scan(file = paste0(input.dir.path, index.files[sss], ".index"))
  file.path = paste0(input.dir.path, index.files[sss], ".info")
  
siteSize=2048
treatment='Copper'
strand='both'

## directory name 
alt.name=paste0(treatment,".", siteSize, ".", strand,  ".alt")
null.name=paste0(treatment,".", siteSize, ".", strand, ".null")
all.name = paste0(treatment,".", siteSize, ".", strand)


## DESeq2 300
deseq.dat = read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.all.pval.0.txt"))

write.table(deseq.dat[index.here,], file= paste0(file.path, ".300"),  col.names = FALSE, row.names = FALSE, quote=FALSE)

## DESeq2 600
deseq.dat = read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".alt.run/output/DESeq2.600.all.pval.0.txt"))

write.table(deseq.dat[index.here,], file= paste0(file.path, ".600"),  col.names = FALSE, row.names = FALSE, quote=FALSE)




## read chromosome and sites number    
path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/info/", treatment, ".", siteSize, ".chr.sites.Robj")
load(file=path)
# chr.list and sites.list
chr.all = unlist(chr.list)
sites.all = unlist(sites.list)


## read logLR from wavelet
## alt
load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/wave/", alt.name, ".run/sum/", alt.name, ".Robj"))
wave.alt = unlist(logLR_list)
wave.done.alt = unlist(done_list)

## read logLR from multiseq
## alt
load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/multiscale/", alt.name, ".sum/", alt.name, ".Robj"))
ms.alt = unlist(logLR_list)
ms.done.alt = unlist(done_list)

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



############################
## collect information
############################

##index.here
chr.here = chr.all[index.here]
sites.here = sites.all[index.here]


### get location information from chr.here and sites.here
st.posi.here = en.posi.here = rep(NA, length(index.here))
for(chr in 1:22){
    wh = which(chr.here == chr)
    path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/", treatment, ".", siteSize, ".", "chr", chr, ".locus")
    locus.info = read.table(path, header = TRUE, as.is = TRUE)
    st.posi.here[wh] = locus.info$st.posi[sites.here[wh]]
    en.posi.here[wh] = locus.info$en.posi[sites.here[wh]]
}    


pval.wave.here = pval.wave[index.here]
pval.ms.here = pval.ms[index.here]
pval.deseq.full.here = pval.deseq.full.0[index.here]
pval.deseq.600.here = pval.deseq.600.0[index.here]
pval.deseq.300.here = pval.deseq.3.0[index.here]
pval.deseq.100.here = pval.deseq.100.0[index.here]

qval.wave.here = qval.wave[index.here]
qval.ms.here = qval.ms[index.here]
qval.deseq.full.here = qval.deseq.full[index.here]
qval.deseq.600.here = qval.deseq.600[index.here]
qval.deseq.300.here = qval.deseq.300[index.here]
qval.deseq.100.here = qval.deseq.100[index.here]


logLR.wave.here = wave.alt[index.here]
logLR.ms.here = ms.alt[index.here]



##sum(all.info$chr != chr.here)
##sum(all.info$sites != sites.here)
##sum(all.info$st.posi != st.posi.here)
##sum(all.info$en.posi != en.posi.here)
##sum(round(all.info$pval.wave,4) != round(pval.wave.here,4))
##sum(round(all.info$pval.ms,4) != round(pval.ms.here,4))
##sum(round(all.info$pval.deseq,4) != round(pval.deseq.300.here,4))
##sum(round(all.info$qval.wave,4) != round(qval.wave.here,4))
##sum(round(all.info$qval.ms,4) != round(qval.ms.here,4))
##sum(round(all.info$qval.deseq,4) != round(qval.deseq.300.here,4))
##sum(round(all.info$logLR.wave,4) != round(logLR.wave.here,4))
##sum(round(all.info$logLR.ms,4) != round(logLR.ms.here,4))







data.info = data.frame(index = index.here,
    chr = chr.here,
    sites = sites.here,
    st.posi =  st.posi.here,
    en.posi = en.posi.here,
    pval.wave = pval.wave.here,
    pval.ms = pval.ms.here,
    pval.deseq.full = pval.deseq.full.here,
    pval.deseq.600 = pval.deseq.600.here,
    pval.deseq.300 = pval.deseq.300.here,
    pval.deseq.100 = pval.deseq.100.here,
    qval.wave = qval.wave.here,
    qval.ms = qval.ms.here,
    qval.deseq.full = qval.deseq.full.here,
    qval.deseq.600 = qval.deseq.600.here,
    qval.deseq.300 = qval.deseq.300.here,
    qval.deseq.100 = qval.deseq.100.here,
    logLR.wave = logLR.wave.here,
    logLR.ms = logLR.ms.here)



write.table(data.info, file= file.path, col.names = TRUE, row.names = FALSE, quote=FALSE)

  
}




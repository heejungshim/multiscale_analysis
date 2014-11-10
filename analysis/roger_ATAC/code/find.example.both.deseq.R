#!/usr/bin/env Rscript

## Aim : This file contains Rscript to find examples of difference multiseq detects, but deseq missed.
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


## read chromosome and sites number    
path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/info/", treatment, ".", siteSize, ".chr.sites.Robj")
load(file=path)
# chr.list and sites.list


## Here copy from see.statistic.3method.both.R
## read p-value from DESeq
#deseq.100.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".alt.run/output/min.pval.txt"))[,1])

#deseq.100.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 100, ".null.run/output/min.pval.txt"))[,1])

deseq.300.alt = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.min.pval.0.txt"))[,1])

deseq.300.null = as.numeric(read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".null.run/output/DESeq2.min.pval.0.txt"))[,1])

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



## read empirical p-values
out.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.pval.discrete.Robj")
load(out.path)
## pval.deseq.3.0

input.path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/pval.ms.wave.", all.name, ".Robj")
load(input.path)
## pval.ms pval.wave






## find tests with invalid outputs and remove them from the analysis
##deseq.100.alt.inf = which(is.infinite(deseq.100.alt)==TRUE)
##deseq.300.alt.inf = which(is.infinite(deseq.300.alt)==TRUE)
##deseq.100.null.inf = which(is.infinite(deseq.100.null)==TRUE)
##deseq.300.null.inf = which(is.infinite(deseq.300.null)==TRUE)

##del.ix = union(union(union(deseq.100.alt.inf, deseq.100.null.inf), which(is.na(ms.alt)==TRUE)), which(is.na(ms.null)==TRUE))

del.ix = union(union(which(is.na(deseq.300.alt) == TRUE), which(is.na(deseq.300.null)==TRUE)), union(which(is.na(pval.ms)==TRUE), which(is.na(pval.wave)== TRUE)))



wave.a = wave.alt[-del.ix]
wave.n = wave.null[-del.ix]
ms.a = ms.alt[-del.ix]
ms.n = ms.null[-del.ix]
##deseq.100.a = deseq.100.alt[-del.ix]
##deseq.300.a = deseq.300.alt[-del.ix]
##deseq.100.n = deseq.100.null[-del.ix]
##deseq.300.n = deseq.300.null[-del.ix]
deseq.300.a = deseq.300.alt[-del.ix]
deseq.300.n = deseq.300.null[-del.ix]




## delete them from chr and sites
chr.all = unlist(chr.list)[-del.ix]
sites.all = unlist(sites.list)[-del.ix]

## get p-value
##load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/deseq.", all.name, ".Robj"))
# pval.deseq.100 and pval.deseq.300
##load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/wave.new.", all.name, ".Robj"))
# pval.wave.new
##load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/summary/ms.new.", all.name, ".Robj"))
# pval.ms.new


p.deseq = pval.deseq.3.0[-del.ix]
p.ms = pval.ms[-del.ix]
p.wave = pval.wave[-del.ix]
index.all = (1:length(pval.ms))[-del.ix]


length(chr.all)
length(sites.all)
length(p.deseq)
length(p.ms)
length(p.wave)
length(index.all)

## [1] 186399
## [1] 186399
## [1] 186399
## [1] 186399
## [1] 186399


#####################
### Storey FDR 
#####################

#install.packages("qvalue")
library("qvalue")

qval.wave = qvalue(p.wave)
qval.ms = qvalue(p.ms)
qval.deseq = qvalue(p.deseq)


#############################
# let's find examples!!!!!
#############################


sum((qval.ms$qvalues < 0.01) & (qval.deseq$qvalues > 0.5))
## 61
sum((qval.ms$qvalues < 0.01) & (qval.deseq$qvalues < 0.02))
## 48


msOnly.index = which((qval.ms$qvalues < 0.01) & (qval.deseq$qvalues > 0.5))
both.index = which((qval.ms$qvalues < 0.01) & (qval.deseq$qvalues < 0.02))

# chr sites st.posi en.posi pval.wave pval.ms qval.wave qval.ms logLR.wave logLR.ms logLR.wave.null logLR.ms.null



length(chr.all)
length(sites.all)
length(p.wave)
length(p.ms)
length(p.deseq)
length(qval.wave$qvalues)
length(qval.ms$qvalues)
length(qval.deseq$qvalues)
length(wave.a)
length(ms.a)
length(wave.n)
length(ms.n)
length(index.all)

sel.ix = msOnly.index

chr.here = chr.all[sel.ix]
sites.here = sites.all[sel.ix]
pval.wave.here = p.wave[sel.ix]
pval.ms.here = p.ms[sel.ix]
pval.deseq.here = p.deseq[sel.ix]
qval.wave.here = qval.wave$qvalues[sel.ix]
qval.ms.here = qval.ms$qvalues[sel.ix]
qval.deseq.here = qval.deseq$qvalues[sel.ix]
logLR.wave.here = wave.a[sel.ix]
logLR.ms.here = ms.a[sel.ix]
logLR.wave.null.here = wave.n[sel.ix]
logLR.ms.null.here = ms.n[sel.ix]
index.here = index.all[sel.ix]


### get location information from chr.here and sites.here
st.posi.here = en.posi.here = rep(NA, length(sel.ix))
for(chr in 1:22){
    wh = which(chr.here == chr)
    path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/", treatment, ".", siteSize, ".", "chr", chr, ".locus")
    locus.info = read.table(path, header = TRUE, as.is = TRUE)
    st.posi.here[wh] = locus.info$st.posi[sites.here[wh]]
    en.posi.here[wh] = locus.info$en.posi[sites.here[wh]]
}    



data.info = data.frame(chr = chr.here, sites = sites.here, st.posi =  st.posi.here,  en.posi = en.posi.here, pval.wave = pval.wave.here, pval.ms = pval.ms.here, pval.deseq = pval.deseq.here, qval.wave = qval.wave.here, qval.ms = qval.ms.here, qval.deseq = qval.deseq.here,  logLR.wave = logLR.wave.here, logLR.ms = logLR.ms.here, logLR.wave.null = logLR.wave.null.here, logLR.ms.null = logLR.ms.null.here, index = index.here)


msOnly.info = data.info

write.table(msOnly.info, file= paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/tmp/",  all.name, ".msOnly.DESeq.info"), col.names = TRUE, row.names = FALSE, quote=FALSE)

deseq.dat = read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.all.pval.0.txt"))

write.table(deseq.dat[msOnly.info$index,], file= paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/tmp/",  all.name, ".msOnly.DESeq.all.pval.info"), col.names = TRUE, row.names = FALSE, quote=FALSE)







sel.ix = both.index


chr.here = chr.all[sel.ix]
sites.here = sites.all[sel.ix]
pval.wave.here = p.wave[sel.ix]
pval.ms.here = p.ms[sel.ix]
pval.deseq.here = p.deseq[sel.ix]
qval.wave.here = qval.wave$qvalues[sel.ix]
qval.ms.here = qval.ms$qvalues[sel.ix]
qval.deseq.here = qval.deseq$qvalues[sel.ix]
logLR.wave.here = wave.a[sel.ix]
logLR.ms.here = ms.a[sel.ix]
logLR.wave.null.here = wave.n[sel.ix]
logLR.ms.null.here = ms.n[sel.ix]
index.here = index.all[sel.ix]


### get location information from chr.here and sites.here
st.posi.here = en.posi.here = rep(NA, length(sel.ix))
for(chr in 1:22){
    wh = which(chr.here == chr)
    path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/locus/", treatment, ".", siteSize, ".", "chr", chr, ".locus")
    locus.info = read.table(path, header = TRUE, as.is = TRUE)
    st.posi.here[wh] = locus.info$st.posi[sites.here[wh]]
    en.posi.here[wh] = locus.info$en.posi[sites.here[wh]]
}    



data.info = data.frame(chr = chr.here, sites = sites.here, st.posi =  st.posi.here,  en.posi = en.posi.here, pval.wave = pval.wave.here, pval.ms = pval.ms.here, pval.deseq = pval.deseq.here, qval.wave = qval.wave.here, qval.ms = qval.ms.here, qval.deseq = qval.deseq.here,  logLR.wave = logLR.wave.here, logLR.ms = logLR.ms.here, logLR.wave.null = logLR.wave.null.here, logLR.ms.null = logLR.ms.null.here, index = index.here)


both.info = data.info

write.table(both.info, file= paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/tmp/",  all.name, ".both.DESeq.info"), col.names = TRUE, row.names = FALSE, quote=FALSE)

deseq.dat = read.table(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/", treatment, ".", siteSize, ".", strand, ".", 300, ".alt.run/output/DESeq2.all.pval.0.txt"))

write.table(deseq.dat[both.info$index,], file= paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/tmp/",  all.name, ".both.DESeq.all.pval.info"), col.names = FALSE, row.names = FALSE, quote=FALSE)





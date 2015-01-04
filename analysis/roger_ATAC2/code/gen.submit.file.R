#!/usr/bin/env Rscript

## Aim : This file contains scripts/funtions to generate submit files for ATAC-seq analysis
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+

# submit all files together
#[shim@spudhead Copper.2048.both.alt]$ for i in *sh; do echo qsub -l h_vmem=1g $i; done
#qsub -l h_vmem=1g ms.1.sh
#qsub -l h_vmem=1g ms.2.sh
#[shim@spudhead Copper.2048.both.alt]$ for i in *sh; do qsub -l h_vmem=1g $i; done
#Your job-array 9385303.1-20:1 ("ms.1.sh") has been submitted
#Your job-array 9385304.1-20:1 ("ms.2.sh") has been submitted




###########################################
# preprocess data for wavelets and DESeq and run multiseq together 
###########################################


com.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/gen.data/com/'
wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/'
siteSize=1024
treatment='Copper'
null=FALSE
strand='both'
meanR.thresh=2
window.size.list=c(100,300,1024)
wavelet.preprocess.QT=TRUE
wavelet.preprocess.NoQT=TRUE
deseq.preprocess=TRUE
get.com.preprocess.run.multiseq.ATACseq(com.path, wd.path, siteSize, treatment, null, strand, meanR.thresh, window.size.list, wavelet.preprocess.QT, wavelet.preprocess.NoQT, deseq.preprocess)



get.com.preprocess.run.multiseq.ATACseq <- function(com.path, wd.path, siteSize, treatment, null, strand, meanR.thresh, window.size.list, wavelet.preprocess.QT, wavelet.preprocess.NoQT, deseq.preprocess){
   
    ## directory name 
    if(!null){
        com.dir.name=paste0(treatment,".", siteSize, ".", strand,  ".alt")
    }else{
        com.dir.name=paste0(treatment,".", siteSize, ".", strand, ".null")
    }
    ## make directory 
    com.out.dir.path = paste0(com.path, com.dir.name) 
    if(!file.exists(com.out.dir.path)){
        dir.create(com.out.dir.path)
    }
    ## make err directory 
    if(!file.exists(paste0(com.out.dir.path, "/err"))){
        dir.create(paste0(com.out.dir.path, "/err"))
    }
    ## read number of sites for each chromosome
    path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/locus/", treatment, ".", siteSize, ".numSites.txt")
    numSites.list = scan(path)

    ## for each chromosome
    for(chr in 1:22){

        file.name = paste0(com.out.dir.path, "/all.", chr, ".sh")

        com = "#!/bin/bash"
        cat(com, file = file.name)
        cat("\n", file = file.name, append = TRUE)

        com = paste("#$ -t 1-", numSites.list[chr], sep="")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)	

        com = paste("#$ -o ", com.out.dir.path, "/err/out.", chr, sep="")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)	

        com = paste("#$ -e ", com.out.dir.path, "/err/err.", chr, sep="")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)	

        com = paste0("cd ", wd.path)
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)

        com = paste0("/data/tools/R-3.1.0/bin/R CMD BATCH --no-save --no-restore \"--args chr=", chr, " sites.ix=$SGE_TASK_ID wd.path='", wd.path, "' siteSize=", siteSize, " treatment='", treatment, "' null=", null, " strand='", strand, "' meanR.thresh=", meanR.thresh, " window.size.list=c(")
        if(length(window.size.list) > 1){
            for(w in 1:(length(window.size.list) -1)){
                com = paste0(com, window.size.list[w], ",")
            }
        }
        com = paste0(com, window.size.list[length(window.size.list)])
        com = paste0(com, ") wavelet.preprocess.QT=", wavelet.preprocess.QT, " wavelet.preprocess.NoQT=", wavelet.preprocess.NoQT, " deseq.preprocess=", deseq.preprocess, "\" /mnt/lustre/home/shim/multiscale_analysis/src/R/preprocess.wave.DESeq.N.run.multiscale.on.roger.ATACseq.R")
        cat(com, file = file.name, append = TRUE)
	cat("\n", file = file.name, append = TRUE)
    }
}






###########################
## collect logLR from multiseq analysis
###########################


collect.logLR.multiseq <- function(com.path, wd.path, siteSize, treatment, null, strand){
    
    ## directory name 
    if(!null){
        com.dir.name=paste0(treatment,".", siteSize, ".", strand, ".alt")
    }else{
        com.dir.name=paste0(treatment,".", siteSize, ".", strand, ".null")
    }
    ## make directory 
    com.out.dir.path = paste0(com.path, com.dir.name, ".post")
    if(!file.exists(com.path)){
        dir.create(com.path)
    }
    if(!file.exists(com.out.dir.path)){
        dir.create(com.out.dir.path)
    }
    ## make err directory 
    if(!file.exists(paste0(com.out.dir.path, "/err"))){
        dir.create(paste0(com.out.dir.path, "/err"))
    }
    
    ## read number of sites for each chromosome
    path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/locus/", treatment, ".", siteSize, ".numSites.txt")
    numSites.list = scan(path)

    for(chr in 1:22){

        numSites = numSites.list[chr]
        file.name = paste0(com.out.dir.path, "/P.ms.", chr, ".sh")

        com = "#!/bin/bash"
        cat(com, file = file.name)
        cat("\n", file = file.name, append = TRUE)

        com = paste("#$ -o ", com.out.dir.path, "/err/out.", chr, ".txt", sep="")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)	

        com = paste("#$ -e ", com.out.dir.path, "/err/err.", chr, ".txt", sep="")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)	

        com = paste0("cd ", wd.path)
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)

        com = paste0("/data/tools/R-3.1.0/bin/R CMD BATCH --no-save --no-restore \"--args chr=", chr, " st.sites=1 en.sites=", numSites, " path.output.dir='", wd.path, com.dir.name, ".output/' path.sum.dir='", wd.path,  com.dir.name, ".sum/' output.file.name='", com.dir.name, "'\" /mnt/lustre/home/shim/multiscale_analysis/src/R/get.logLR.ms.R")    
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)
    }
}



com.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/multiscale/com/'
wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/multiscale/'
siteSize=1024
treatment='Copper'
null=FALSE
strand='both'


##############################################
# Check if we have a complete list of logLR
##############################################
siteSize=1024
treatment='Copper'
null=FALSE
strand='both'
wd.path = '/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/multiscale/'
if(!null){
  com.dir.name=paste0(treatment,".", siteSize, ".", strand, ".alt")
}else{
  com.dir.name=paste0(treatment,".", siteSize, ".", strand, ".null")
}
output.path = paste0(wd.path, com.dir.name, ".sum/", com.dir.name)
numF = rep(NA, 22)
for(chr in 1:22){
  load(paste0(output.path, ".", chr, ".Robj"))
  numF[chr] = sum(!done)
}
numF






###########################
## run wavelets with QT
###########################
com.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/wave/com/'
wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/wave/'
siteSize=1024
treatment='Copper'
null=FALSE
strand='both'

run.wavelets <- function(com.path, wd.path, siteSize, treatment, null, strand){
    
    ## directory name 
    if(!null){
        com.dir.name=paste0(treatment,".", siteSize, ".", strand, ".alt")
    }else{
        com.dir.name=paste0(treatment,".", siteSize, ".", strand, ".null")
    }
    ## make directory
    if(!file.exists(com.path)){
      dir.create(com.path)
    }
    com.out.dir.path = paste0(com.path, com.dir.name) 
    if(!file.exists(com.out.dir.path)){
        dir.create(com.out.dir.path)
    }
    ## make err directory 
    if(!file.exists(paste0(com.out.dir.path, "/err"))){
        dir.create(paste0(com.out.dir.path, "/err"))
    }

    ## make output directory 
    wave.out.dir.path = paste0(wd.path, com.dir.name, ".run") 
    if(!file.exists(wave.out.dir.path)){
        dir.create(wave.out.dir.path)
    }

    ## read number of sites for each chromosome
    path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/locus/", treatment, ".", siteSize, ".numSites.txt")
    numSites.list = scan(path)

    for(chr in 1:22){

        numSites = numSites.list[chr]
        file.name = paste0(com.out.dir.path, "/wave.", chr, ".sh")

        com = "#!/bin/bash"
        cat(com, file = file.name)
        cat("\n", file = file.name, append = TRUE)

        com = paste("#$ -t 1-", numSites, sep="")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)	

        com = paste("#$ -o ", com.out.dir.path, "/err/out.", chr, ".txt", sep="")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)	

        com = paste("#$ -e ", com.out.dir.path, "/err/err.", chr, ".txt", sep="")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)	

        com = paste0("cd ", wd.path, com.dir.name, ".run/")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)

        com = paste0("~shim/WaveQTL/WaveQTL -gmode 1 -g /mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/info/geno6.txt -p /mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/wave/", com.dir.name, ".data/WC.", chr, ".$SGE_TASK_ID.txt -u /mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/wave/", com.dir.name, ".data/use.", chr, ".$SGE_TASK_ID.txt -o ", com.dir.name, ".", chr, ".$SGE_TASK_ID -f ", siteSize," -fph 1 -nullcheck 1")   
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)
    }
}


###########################
## run wavelets without  QT
###########################
com.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/waveNoQT/com/'
wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/waveNoQT/'
siteSize=1024
treatment='Copper'
null=FALSE
strand='both'
run.wavelets.noqt <- function(com.path, wd.path, siteSize, treatment, null, strand){
    
    ## directory name 
    if(!null){
        com.dir.name=paste0(treatment,".", siteSize, ".", strand, ".alt")
    }else{
        com.dir.name=paste0(treatment,".", siteSize, ".", strand, ".null")
    }
    ## make directory
    if(!file.exists(com.path)){
      dir.create(com.path)
    }
    com.out.dir.path = paste0(com.path, com.dir.name) 
    if(!file.exists(com.out.dir.path)){
        dir.create(com.out.dir.path)
    }
    ## make err directory 
    if(!file.exists(paste0(com.out.dir.path, "/err"))){
        dir.create(paste0(com.out.dir.path, "/err"))
    }

    ## make output directory 
    wave.out.dir.path = paste0(wd.path, com.dir.name, ".run") 
    if(!file.exists(wave.out.dir.path)){
        dir.create(wave.out.dir.path)
    }

    ## read number of sites for each chromosome
    path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/locus/", treatment, ".", siteSize, ".numSites.txt")
    numSites.list = scan(path)

    for(chr in 1:22){

        numSites = numSites.list[chr]
        file.name = paste0(com.out.dir.path, "/waveNoQT.", chr, ".sh")

        com = "#!/bin/bash"
        cat(com, file = file.name)
        cat("\n", file = file.name, append = TRUE)

        com = paste("#$ -t 1-", numSites, sep="")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)	

        com = paste("#$ -o ", com.out.dir.path, "/err/out.", chr, ".txt", sep="")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)	

        com = paste("#$ -e ", com.out.dir.path, "/err/err.", chr, ".txt", sep="")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)	

        com = paste0("cd ", wd.path, com.dir.name, ".run/")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)

        com = paste0("~shim/WaveQTL/WaveQTL -gmode 1 -g /mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/info/geno6.txt -p /mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/waveNoQT/", com.dir.name, ".data/WC.", chr, ".$SGE_TASK_ID.txt -u /mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/waveNoQT/", com.dir.name, ".data/use.", chr, ".$SGE_TASK_ID.txt -o ", com.dir.name, ".", chr, ".$SGE_TASK_ID -f ", siteSize," -fph 1 -nullcheck 1")   
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)
    }
}



###########################
## collect logLR from wavelet analysis
###########################
## with QT
com.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/wave/com/'
wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/wave/'
siteSize=1024
treatment='Copper'
null=FALSE
strand='both'
## without QT
com.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/waveNoQT/com/'
wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/waveNoQT/'
siteSize=1024
treatment='Copper'
null=FALSE
strand='both'

collect.logLR.wavelets <- function(com.path, wd.path, siteSize, treatment, null, strand){
    
    ## directory name 
    if(!null){
        com.dir.name=paste0(treatment,".", siteSize, ".", strand, ".alt")
    }else{
        com.dir.name=paste0(treatment,".", siteSize, ".", strand, ".null")
    }
    ## make directory 
    com.out.dir.path = paste0(com.path, com.dir.name, ".post") 
    if(!file.exists(com.out.dir.path)){
        dir.create(com.out.dir.path)
    }
    ## make err directory 
    if(!file.exists(paste0(com.out.dir.path, "/err"))){
        dir.create(paste0(com.out.dir.path, "/err"))
    }
    
    ## read number of sites for each chromosome
    path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/locus/", treatment, ".", siteSize, ".numSites.txt")
    numSites.list = scan(path)

    for(chr in 1:22){

        numSites = numSites.list[chr]
        file.name = paste0(com.out.dir.path, "/P.wave.", chr, ".sh")

        com = "#!/bin/bash"
        cat(com, file = file.name)
        cat("\n", file = file.name, append = TRUE)

        com = paste("#$ -o ", com.out.dir.path, "/err/out.", chr, ".txt", sep="")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)	

        com = paste("#$ -e ", com.out.dir.path, "/err/err.", chr, ".txt", sep="")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)	

        com = paste0("cd ", wd.path, com.dir.name, ".run/")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)

        com = paste0("/data/tools/R-3.1.0/bin/R CMD BATCH --no-save --no-restore \"--args chr=", chr, " st.sites=1 en.sites=", numSites, " path.output.dir='", wd.path, com.dir.name, ".run/' output.file.name='", com.dir.name, "'\" /mnt/lustre/home/shim/multiscale_analysis/src/R/get.logLR.wave.R")    
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)
    }
}







###########################
## Collect data for DESeq
###########################



com.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/deseq/com/'
wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/deseq/'
siteSize=1024
treatment='Copper'
null=FALSE
strand='both'
window.size=100
#window.size=300
#window.size=1024
numSam=6

combine.for.DESeq.data <- function(com.path, wd.path, siteSize, treatment, null, strand, window.size, numSam){
    
    
    ## directory name 
    if(!null){
        com.dir.name=paste0(treatment,".", siteSize, ".", strand,  ".", window.size, ".alt")
    }else{
        com.dir.name=paste0(treatment,".", siteSize, ".", strand, ".", window.size, ".null")
    }
    ## make directory 
    com.out.dir.path = paste0(com.path, com.dir.name) 
    if(!file.exists(com.out.dir.path)){
        dir.create(com.out.dir.path)
    }
    ## make err directory 
    if(!file.exists(paste0(com.out.dir.path, "/err"))){
        dir.create(paste0(com.out.dir.path, "/err"))
    }
    ## read number of sites for each chromosome
    path = paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/locus/", treatment, ".", siteSize, ".numSites.txt")
    numSites.list = scan(path)

    for(chr in 1:22){

        numSites = numSites.list[chr]
        file.name = paste0(com.out.dir.path, "/deseq.", chr, ".sh")

        com = "#!/bin/bash"
        cat(com, file = file.name)
        cat("\n", file = file.name, append = TRUE)

        com = paste("#$ -o ", com.out.dir.path, "/err/out.", chr, ".txt", sep="")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)	

        com = paste("#$ -e ", com.out.dir.path, "/err/err.", chr, ".txt", sep="")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)	

        com = paste0("cd ", wd.path)
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)

        com = paste0("/data/tools/R-3.1.0/bin/R CMD BATCH --no-save --no-restore \"--args chr=", chr, " st.sites=1 en.sites=", numSites, " wd.path='", wd.path, "' siteSize=", siteSize, " treatment='", treatment, "' null=", null, " strand='", strand, "' window.size=", window.size, " numSam=6\" /mnt/lustre/home/shim/multiscale_analysis/src/R/combine.for.DESeq.data.R")
        cat(com, file = file.name, append = TRUE)
        cat("\n", file = file.name, append = TRUE)
    }
}








###########################
## run DESeq
###########################

wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/deseq/'
siteSize=1024
treatment='Copper'
null=FALSE
strand='both'
window.size=100
#window.size=300
#window.size=1024
numSam=6
filter.cut=0

run.DESeq <- function(wd.path, siteSize, treatment, null, strand, window.size, numSam){
        
    ## directory name 
    if(!null){
        com.dir.name=paste0(treatment,".", siteSize, ".", strand,  ".", window.size, ".alt")
    }else{
        com.dir.name=paste0(treatment,".", siteSize, ".", strand, ".", window.size, ".null")
    }

    ## make directory 
    com.out.dir.path = paste0(wd.path, com.dir.name, ".run/com") 
    if(!file.exists(com.out.dir.path)){
        dir.create(com.out.dir.path)
    }
    ## make err directory 
    if(!file.exists(paste0(com.out.dir.path, "/err"))){
        dir.create(paste0(com.out.dir.path, "/err"))
    }

    file.name = paste0(com.out.dir.path, "/DESeq.sh")

    com = "#!/bin/bash"
    cat(com, file = file.name)
    cat("\n", file = file.name, append = TRUE)
    
    com = paste("#$ -o ", com.out.dir.path, "/err/out.txt", sep="")
    cat(com, file = file.name, append = TRUE)
    cat("\n", file = file.name, append = TRUE)	

    com = paste("#$ -e ", com.out.dir.path, "/err/err.txt", sep="")
    cat(com, file = file.name, append = TRUE)
    cat("\n", file = file.name, append = TRUE)	

    com = paste0("cd ", com.out.dir.path)
    cat(com, file = file.name, append = TRUE)
    cat("\n", file = file.name, append = TRUE)

    com = paste0("/data/tools/R-3.1.0/bin/R CMD BATCH --no-save --no-restore \"--args wd.path='", wd.path, "' siteSize=", siteSize, " treatment='", treatment, "' null=", null, " strand='", strand, "' window.size=", window.size, " numSam=", numSam, " filter.cut=",filter.cut, "\" /mnt/lustre/home/shim/multiscale_analysis/src/R/run.DESeq2.on.roger.ATACseq.R")    
    cat(com, file = file.name, append = TRUE)
    cat("\n", file = file.name, append = TRUE)

}














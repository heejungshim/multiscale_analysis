#!/usr/bin/env Rscript

## Aim : This file contains Rscript to estimate effect size using multiseq. I modified a script `~/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v2/prepareData/script/est.ES.multiseq.for.simulation.R'.
##
##
## Usage R CMD BATCH --no-save --no-restore "--args ss=$SGE_TASK_ID" est.ES.multiseq.R
## See /mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_simple_v3/com/ES.multiseq.sh
##
## 
## Copyright (C) 2015 Heejung Shim
##
## License: GPL3+

## ss = 1
args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))

setwd("/mnt/lustre/home/shim/multiscale_analysis")
multiscale.analysis.repodir <- scan(".multiscale_analysis.repodir.txt", what=character())
source(paste0(multiscale.analysis.repodir, "/src/R/my.utils.R"))
library("multiseq")
library("ashr")

## set working directory
wd.path = paste0(multiscale.analysis.repodir, "/analysis/simulation/sample_size/simulation_simple_v3/prepareData/")
setwd(wd.path)

## path to output
path.output = "~/multiscale_analysis/analysis/simulation/sample_size/simulation_simple_v3/data/"

## Read phenotype and genotype data
path.data = "~/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v1/data/"
phenoD = as.matrix(read.table(paste0(path.data, "pheno.dat.", ss), as.is = TRUE))
genoD = scan(paste0(path.data, "orig.geno.dat.", ss), what=double())
genoR = as.numeric(round(genoD))

wh0 = which(genoR == 0)
wh1 = which(genoR == 1)
wh2 = which(genoR == 2)
##length(wh0)
##length(wh1)
##length(wh2)


## estimate effect size when we have enough samples with genotypes 0 and 1 
if((length(wh0) > 2) & (length(wh1) > 2)){

  numPerm=NULL
  numSig=10

  ## Use only two genotype classes
  if(length(wh2) > 0){
    genoF = genoR[-wh2]
    phenoF = phenoD[-wh2,]
  }else{
    genoF = genoR
    phenoF = phenoD
  }
    
  ## run multiseq
  res.ES = multiseq(x = as.matrix(phenoF), g = genoF, read.depth = NULL)
  if(length(res.ES$effect.mean) == 1024){

    sig0.smooth = exp(res.ES$baseline.mean)
    
    multiseq.mean = res.ES$effect.mean
    multiseq.sd = sqrt(res.ES$effect.var)
   
    wh = which(abs(multiseq.mean) < 3*multiseq.sd)
    multiseq.mean.sig3 = multiseq.mean
    multiseq.mean.sig3[wh] = 0

    sig1.smooth = exp(res.ES$baseline.mean + multiseq.mean.sig3)

    ## for null signal
    sig.mat = cbind(sig0.smooth, sig1.smooth)
    null.sig = apply(sig.mat, 1, mean)
  
    ## print out signals
    ## alt
    write.table(sig0.smooth, file = paste0(path.output, "alt.sig0.", ss), col.names = FALSE, row.names = FALSE, quote=FALSE)
    write.table(sig1.smooth, file = paste0(path.output, "alt.sig1.", ss), col.names = FALSE, row.names = FALSE, quote=FALSE)
    ## null
    write.table(null.sig, file = paste0(path.output, "null.sig0.", ss), col.names = FALSE, row.names = FALSE, quote=FALSE)
    write.table(null.sig, file = paste0(path.output, "null.sig1.", ss), col.names = FALSE, row.names = FALSE, quote=FALSE)
  }else{
      cat("", file = paste0(path.output, "NG.", ss))
    }

}else{
  cat("", file = paste0(path.fig, "NG.", ss))
}






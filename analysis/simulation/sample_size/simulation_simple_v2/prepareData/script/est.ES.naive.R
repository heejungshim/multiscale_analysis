#!/usr/bin/env Rscript

## Aim : This file contains Rscript to estimate effect size using BAYES.THR in wavethresh package. I modified a script `~/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v1/prepareData/script/est.ES.naive.for.simulation.R'.
##
##
## Usage R CMD BATCH --no-save --no-restore "--args ss=$SGE_TASK_ID" est.ES.naive.R
## See /mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_simple_v2/com/ES.naive.sh
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

## set working directory
wd.path = paste0(multiscale.analysis.repodir, "/analysis/simulation/sample_size/simulation_simple_v2/prepareData/")
setwd(wd.path)
## output path
path.output = "~/multiscale_analysis/analysis/simulation/sample_size/simulation_simple_v2/data/"


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

  ## Take a mean profile
  sig0 = apply(phenoD[wh0,], 2, mean)
  sig1 = apply(phenoD[wh1,], 2, mean)

  ## denoise mean curves using wavethresh
  library(wavethresh)

  ## combine two signals
  sig.all = c(sig0, sig1)
  sig.all.smooth = BAYES.THR(sig.all)

  sig0.smooth = sig.all.smooth[1:1024]
  sig1.smooth = sig.all.smooth[1025:2048]

  val = 6
  cut.thresh = val/70
  delix = which(abs(sig0.smooth - sig1.smooth) <= cut.thresh)
  sig1.smooth[delix] = sig0.smooth[delix]

  ## handle 0 or negative count
  wh.zero = which(sig0.smooth <= 1/70)
  if(length(wh.zero) > 0){ 
    sig0.smooth[wh.zero] = 0
  }
  wh.zero = which(sig1.smooth <=  1/70)
  if(length(wh.zero) > 0){ 
    sig1.smooth[wh.zero] = 0
  }

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




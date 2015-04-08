#!/usr/bin/env Rscript

## Aim : This file contains Rscript to estimate effect size using BAYES.THR in wavethresh package. I modified a script `~/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/code/try.different.effect.size.R'.
##
##
## Usage R CMD BATCH --no-save --no-restore "--args ss=$SGE_TASK_ID" est.ES.naive.for.simulation.R
## See /mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v1/com/ES.naive.sh
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
##WaveQTL.repodir <- scan(".WaveQTL.repodir.txt", what=character())

## set working directory
wd.path = paste0(multiscale.analysis.repodir, "/analysis/simulation/sample_size/simulation_manydsQTL_v1/prepareData/")
setwd(wd.path)

## Read phenotype and genotype data
path.fig = "~/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v1/fig/ES/" 
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
    sig0.smooth[wh.zero] = 1/70
  }
  wh.zero = which(sig1.smooth <=  1/70)
  if(length(wh.zero) > 0){ 
    sig1.smooth[wh.zero] = 1/70
  }
  
  ## get smooth signal as a ratio
  smooth.ratio = sig1.smooth/sig0.smooth
  write.table(smooth.ratio, file = paste0(path.data, "smooth.ratio.", ss), col.names = FALSE, row.names = FALSE, quote=FALSE)

##########################################
## Make a figure to see how they look like
##########################################
  png(paste0(path.fig, "ESnaive", ss ,".png"), height = 4, width = 7, units="in", res=300)

  numBPs = 1024
  xmin = 1
  xmax = numBPs 
  xval = xmin:xmax

  nf <- layout(matrix(1:4,4,1,byrow = TRUE))

  ## Group 0
  ymax = max(sig0, sig0.smooth, sig1, sig1.smooth)
  par(mar=c(3,3,1,1))
  plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(0, ymax),axes= FALSE, main ="Major Homozygotes")
  lines(xval, sig0, type = "l", col = "orange")
  lines(xval, sig0.smooth, type = "l", col = "red")
  axis(2, font=2)
  box()

  ## Group 1
  ymax = max(sig0, sig0.smooth, sig1, sig1.smooth)
  par(mar=c(3,3,1,1))
  plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(0, ymax),axes= FALSE, main ="Heterozygotes")
  lines(xval, sig1, type = "l", col = "skyblue")
  lines(xval, sig1.smooth, type = "l", col = "blue")
  axis(2, font=2)
  box()

  ## difference
  ymax = max(sig0.smooth - sig1.smooth)
  ymin = min(sig0.smooth - sig1.smooth)

  par(mar=c(3,3,1,1))
  plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(ymin, ymax),axes= FALSE, main ="Major Homozygotes - Heterozygotes")
  lines(xval, sig0.smooth - sig1.smooth, type = "l", col = "green")
  abline(h = 0)
  axis(2, font=2)
  box()

  ## ratio 1 
  ymax = max(log(sig0.smooth) - log(sig1.smooth))
  ymin = min(log(sig0.smooth) - log(sig1.smooth))

  par(mar=c(3,3,1,1))
  plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(ymin, ymax),axes= FALSE, main ="log(Major Homozygotes) - log(Heterozygotes)")
  lines(xval, log(sig0.smooth) - log(sig1.smooth), type = "l", col = "green")
  abline(h = 0)
  axis(2, font=2)
  box()

  dev.off()

  png(paste0(path.fig, "ESsimunaive", ss ,".png"), height = 4, width = 7, units="in", res=300)

  nf <- layout(matrix(1:4,4,1,byrow = TRUE))

  ## read raw data from which we will sample
  raw.data = phenoD
  raw.data.T = ceiling(as.numeric(apply(raw.data, 2, sum)))

  mu0.sig = rep(1/70, 1024)
  mu1.sig = mu0.sig*smooth.ratio
    
  ## ok!!
  ## upper and lower bound!
  trunc.fun = function(x){
    x = max(0, x)
    return(min(1,x))
  }
  mu0.sig = sapply(mu0.sig, trunc.fun)
  mu1.sig = sapply(mu1.sig, trunc.fun)

  sig0 = mu0.sig*raw.data.T
  sig1 = mu1.sig*raw.data.T

  ## Data 
  ymax = max(raw.data.T)
  par(mar=c(1,3,1,1))
  plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(0, ymax),axes= FALSE, main ="Raw data")
  lines(xval, raw.data.T, type = "l", col = "black")
  axis(2, font=2)
  box()

  ## Expected mean for two groups
  ymax = max(sig0, sig1)
  par(mar=c(1,3,1,1))
  plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(0, ymax),axes= FALSE, main ="Expected mean g0: red g1: blue [original]")
  lines(xval, sig0, type = "l", col = "red")
  lines(xval, sig1, type = "l", col = "blue")
  axis(2, font=2)
  box()

  ## Expected difference
  ymax = max(sig1 - sig0)
  ymin = min(sig1 - sig0)
  par(mar=c(1,3,1,1))
  plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(ymin, ymax),axes= FALSE, main ="Expected difference (g1 - g0) [original]")
  lines(xval, sig1 - sig0, type = "l", col = "green")
  abline(h = 0)
  axis(2, font=2)
  box()

  ## Expected log ratio
  wh1 = which(sig1 == 0)
  wh0 = which(sig0 == 0)
  min.val = min(sig1[-wh1], sig0[-wh0]) 
  sig1[wh1] = min.val
  sig0[wh0] = min.val

  ymax = max(log(sig1) - log(sig0))
  ymin = min(log(sig1) - log(sig0))
  par(mar=c(1,3,1,1))
  plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(ymin, ymax),axes= FALSE, main ="Expected log ratio (log(g1) - log(g0)) [original]")
  lines(xval, log(sig1) - log(sig0), type = "l", col = "green")
  abline(h = 0)
  axis(2, font=2)
  box()

  dev.off()

}else{
  cat("", file = paste0(path.fig, "NG.", ss))
}




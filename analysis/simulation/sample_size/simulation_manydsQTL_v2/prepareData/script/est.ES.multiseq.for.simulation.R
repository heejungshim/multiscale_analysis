#!/usr/bin/env Rscript

## Aim : This file contains Rscript to estimate effect size using multiseq. I modified a script in the section `Estimating effect size from multiseq or WaveQTL' of the document `~/mydoc/multiscale/simu_WaveQTL_multiseq/ESfrom2Methods.org'.
##
##
## Usage R CMD BATCH --no-save --no-restore "--args ss=$SGE_TASK_ID" est.ES.multiseq.for.simulation.R
## See /mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v2/com/ES.multiseq.sh
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
wd.path = paste0(multiscale.analysis.repodir, "/analysis/simulation/sample_size/simulation_manydsQTL_v2/prepareData/")
setwd(wd.path)

## WaveQTL path
path.WaveQTL = '/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v3/WaveQTL/'

## path to output
path.fig = "~/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v2/fig/ES/"
path.output = "~/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v2/data/"

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

  ## To make a figure later
  ## Take a mean profile
  sig0 = apply(phenoD[wh0,], 2, mean)
  sig1 = apply(phenoD[wh1,], 2, mean)

  ## combine two signals
  library(wavethresh)
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
    
  ## Save genotype and phenotype data for WaveQTL
  write.table(phenoF, file = paste0(path.WaveQTL, "/phenoD.", ss), row.names=FALSE, col.names=FALSE, quote=FALSE)

  path.this = paste0(path.WaveQTL, "/genoD.", ss)
  cat("rs A T ", file = path.this)
  cat(genoF, file = path.this, append = TRUE)

  ## run multiseq
  res.ES = multiseq(x = as.matrix(phenoF), g = genoF, read.depth = NULL)
  if(length(res.ES$effect.mean) == 1024){
    
    multiseq.mean = res.ES$effect.mean
    multiseq.sd = sqrt(res.ES$effect.var)
    multiseq.base.mean = res.ES$baseline.mean

    wh = which(abs(multiseq.mean) < 2*multiseq.sd)
    multiseq.mean.sig2 = multiseq.mean
    multiseq.mean.sig2[wh] = 0

    wh = which(abs(multiseq.mean) < 3*multiseq.sd)
    multiseq.mean.sig3 = multiseq.mean
    multiseq.mean.sig3[wh] = 0

    ## Save effect size estimate
    this.path = paste0(path.output, "smooth.ratio.2.", ss)
    write.table(exp(multiseq.mean.sig2), file = this.path, row.names=FALSE, col.names = FALSE, quote = FALSE)
    this.path = paste0(path.output, "smooth.ratio.3.", ss)
    write.table(exp(multiseq.mean.sig3), file = this.path, row.names=FALSE, col.names = FALSE, quote = FALSE)

    ## Make a figure to see how they look like
    png(paste0(path.fig, "ESmultiseq", ss ,".png"), height = 4, width = 7, units="in", res=300)

    numBPs = 1024
    xmin = 1
    xmax = numBPs 
    xval = xmin:xmax

    nf <- layout(matrix(1:3,3,1,byrow = TRUE))

    ## Group 0
    ymax = max(sig0, sig0.smooth, sig1, sig1.smooth, exp(multiseq.base.mean))
    par(mar=c(3,3,1,1))
    plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(0, ymax),axes= FALSE, main ="Major Homozygotes, dark green from multiseq")
    lines(xval, sig0, type = "l", col = "orange")
    lines(xval, sig0.smooth, type = "l", col = "red")
    lines(xval, exp(multiseq.base.mean), type = "l", col = "darkgreen")
    axis(2, font=2)
    box()

    ## Group 1
    ymax = max(sig0, sig0.smooth, sig1, sig1.smooth, exp(multiseq.base.mean))
    par(mar=c(3,3,1,1))
    plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(0, ymax),axes= FALSE, main ="Heterozygotes")
    lines(xval, sig1, type = "l", col = "skyblue")
    lines(xval, sig1.smooth, type = "l", col = "blue")
    axis(2, font=2)
    box()

    ## ratio 1
    ymax = max(log(sig1.smooth) - log(sig0.smooth), multiseq.mean)
    ymin = min(log(sig1.smooth) - log(sig0.smooth), multiseq.mean)

    par(mar=c(3,3,1,1))
    plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(ymin, ymax),axes= FALSE, main ="log(Heterozygotes) - log(Major Homozygotes), effect size from multiseq (0, 2, 3 sd)")
    lines(xval, log(sig1.smooth) - log(sig0.smooth), type = "l", col = "green")
    lines(xval, multiseq.mean, type = "l", col = "orange")
    lines(xval, multiseq.mean.sig2, type = "l", col = "red")
    lines(xval, multiseq.mean.sig3, type = "l", col = "blue")
    abline(h = 0)
    axis(2, font=2)
    box()

    dev.off()


    png(paste0(path.fig, "ESsimumultiseq", ss ,".png"), height = 4, width = 7, units="in", res=300)

    nf <- layout(matrix(1:4,4,1,byrow = TRUE))

    ## read raw data from which we will sample
    raw.data = phenoD
    raw.data.T = ceiling(as.numeric(apply(raw.data, 2, sum)))

    mu0.sig = rep(1/70, 1024)
    mu1.sig = mu0.sig*exp(multiseq.mean.sig3)
    
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
    plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(0, ymax),axes= FALSE, main ="Expected mean g0: red g1: blue [multiseq]")
    lines(xval, sig0, type = "l", col = "red")
    lines(xval, sig1, type = "l", col = "blue")
    axis(2, font=2)
    box()

    ## Expected difference
    ymax = max(sig1 - sig0)
    ymin = min(sig1 - sig0)
    par(mar=c(1,3,1,1))
    plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(ymin, ymax),axes= FALSE, main ="Expected difference (g1 - g0) [multiseq]")
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
    plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(ymin, ymax),axes= FALSE, main ="Expected log ratio (log(g1) - log(g0)) [multiseq]")
    lines(xval, log(sig1) - log(sig0), type = "l", col = "green")
    abline(h = 0)
    axis(2, font=2)
    box()
    dev.off()

  }else{
    cat("", file = paste0(path.fig, "NG.", ss))
  }

}else{
  cat("", file = paste0(path.fig, "NG.", ss))
}






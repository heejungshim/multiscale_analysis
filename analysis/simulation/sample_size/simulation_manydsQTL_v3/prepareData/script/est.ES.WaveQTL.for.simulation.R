#!/usr/bin/env Rscript

## Aim : This file contains Rscript to estimate effect size using multiseq. I modified a script in the section `Estimating effect size from multiseq or WaveQTL' of the document `~/mydoc/multiscale/simu_WaveQTL_multiseq/ESfrom2Methods.org'.
##
##
## Usage R CMD BATCH --no-save --no-restore "--args ss=$SGE_TASK_ID" est.ES.WaveQTL.for.simulation.R
## See /mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v3/com/ES.WaveQTL.sh
##
## 
## Copyright (C) 2015 Heejung Shim
##
## License: GPL3+

## ss = 1
args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))

setwd("~/WaveQTL/R")
## Read Haar DWT matrix
Wmat_1024 = read.table("../data/DWT/Wmat_1024",as.is = TRUE)
W2mat_1024 = Wmat_1024*Wmat_1024

## Read output from WaveQTL
WaveQTL.path='~/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v3/WaveQTL/'

## set working directory
wd.path = '~/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v3/prepareData/'
setwd(wd.path)

## path to output
path.fig = "~/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v3/fig/ES/"
path.output = "~/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v3/data/"

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

  ## Read output from WaveQTL
  beta_mean_path = paste0(WaveQTL.path, "output/res.", ss, ".fph.mean.txt")
  beta_var_path = paste0(WaveQTL.path, "output/res.", ss, ".fph.var.txt") 
  ## We'll look at effect size of the first SNP in genotype file
  sel_geno_IX = 1
  ## mean 
  beta_mean = as.numeric(read.table(beta_mean_path)[sel_geno_IX,2:1025])
  beta_dataS = as.vector(-matrix(data=beta_mean, nr = 1, nc = 1024)%*%as.matrix(Wmat_1024))
  ## sd
  beta_var = as.numeric(read.table(beta_var_path)[sel_geno_IX,2:1025])
  beta_var_dataS = as.vector(matrix(data=beta_var, nr=1, nc=1024)%*%as.matrix(W2mat_1024))
  beta_sd_dataS = sqrt(beta_var_dataS)

  
  WaveQTL.mean = beta_dataS
  WaveQTL.sd = beta_sd_dataS

  wh = which(abs(WaveQTL.mean) < 2*WaveQTL.sd)
  WaveQTL.mean.sig2 = WaveQTL.mean
  WaveQTL.mean.sig2[wh] = 0

  wh = which(abs(WaveQTL.mean) < 3*WaveQTL.sd)
  WaveQTL.mean.sig3 = WaveQTL.mean
  WaveQTL.mean.sig3[wh] = 0

  raw.data = phenoD
  raw.data.T = ceiling(as.numeric(apply(raw.data, 2, sum)))
  wh = which(raw.data.T == 0)

  this.path = paste0(path.output, "smooth.ratio.2.", ss)
  res = 1 + 70*WaveQTL.mean.sig2/raw.data.T
  res[wh] = 1
  write.table(res, file = this.path, row.names=FALSE, col.names = FALSE, quote = FALSE)

  this.path = paste0(path.output, "smooth.ratio.3.", ss)
  res = 1 + 70*WaveQTL.mean.sig3/raw.data.T
  res[wh] = 1
  write.table(res, file = this.path, row.names=FALSE, col.names = FALSE, quote = FALSE)
  
  ## Make a figure to see how they look like
  png(paste0(path.fig, "ESWaveQTL", ss ,".png"), height = 4, width = 7, units="in", res=300)
  numBPs = 1024
  xmin = 1
  xmax = numBPs 
  xval = xmin:xmax
  
  nf <- layout(matrix(1:3,3,1,byrow = TRUE))
  
  ## Group 0
  ymax = max(sig0, sig0.smooth, sig1, sig1.smooth)
  par(mar=c(3,3,1,1))
  plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(0, ymax),axes= FALSE, main ="Major Homozygotes, dark green from multiseq")
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

  #### difference
  ymax = max(sig1.smooth - sig0.smooth, WaveQTL.mean.sig3, WaveQTL.mean)
  ymin = min(sig1.smooth - sig0.smooth, WaveQTL.mean.sig3, WaveQTL.mean)

  par(mar=c(3,3,1,1))
  plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(ymin, ymax),axes= FALSE, main ="Major Homozygotes - Heterozygotes, effect size from WaveQTL (0, 2, 3 sd)")
  lines(xval, sig1.smooth - sig0.smooth, type = "l", col = "green")
  lines(xval, WaveQTL.mean, type = "l", col = "orange")
  lines(xval, WaveQTL.mean.sig2, type = "l", col = "red")
  lines(xval, WaveQTL.mean.sig3, type = "l", col = "blue")
  abline(h = 0)
  axis(2, font=2)
  box()

  dev.off()

  png(paste0(path.fig, "ESsimuWaveQTL", ss ,".png"), height = 4, width = 7, units="in", res=300)

  nf <- layout(matrix(1:4,4,1,byrow = TRUE))
  
  ## read raw data from which we will sample
  raw.data = phenoD
  raw.data.T = ceiling(as.numeric(apply(raw.data, 2, sum)))

  mu0.sig = rep(1/70, 1024)
  sig0.WaveQTL = mu0.sig*raw.data.T
  sig1.WaveQTL = sig0.WaveQTL - WaveQTL.mean.sig3
  wh = which(sig1.WaveQTL < 0)
  sig1.WaveQTL[wh] = 0

  sig0 = sig0.WaveQTL
  sig1 = sig1.WaveQTL

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
  plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(0, ymax),axes= FALSE, main ="Expected mean g0: red g1: blue [WaveQTL]")
  lines(xval, sig0, type = "l", col = "red")
  lines(xval, sig1, type = "l", col = "blue")
  axis(2, font=2)
  box()

  ## Expected difference
  ymax = max(sig1 - sig0)
  ymin = min(sig1 - sig0)
  par(mar=c(1,3,1,1))
  plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(ymin, ymax),axes= FALSE, main ="Expected difference (g1 - g0) [WaveQTL]")
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
  plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(ymin, ymax),axes= FALSE, main ="Expected log ratio (log(g1) - log(g0)) [WaveQTL]")
  lines(xval, log(sig1) - log(sig0), type = "l", col = "green")
  abline(h = 0)
  axis(2, font=2)
  box()
  dev.off()

}else{
  cat("", file = paste0(path.fig, "NG.", ss))
}

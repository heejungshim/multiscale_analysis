#!/usr/bin/env Rscript

## Aim : This file contains Rscript to preprocess data for WaveQTL analysis to estimate effect size. I modified a script in the section `Estimating effect size from multiseq or WaveQTL' of the document `~/mydoc/multiscale/simu_WaveQTL_multiseq/ESfrom2Methods.org'.
##
##
## Usage R CMD BATCH --no-save --no-restore "--args ss=$SGE_TASK_ID" preprocess.WaveQTL.for.simulation.R
## See /mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v2/com/preprocess.WaveQTL.sh
##
## 
## Copyright (C) 2015 Heejung Shim
##
## License: GPL3+

## ss = 1
args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))

setwd("/mnt/lustre/home/shim/multiscale_analysis")
WaveQTL.repodir <- scan(".WaveQTL.repodir.txt", what=character())
wd.path='~/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v3/WaveQTL'
setwd(wd.path)
source(paste0(WaveQTL.repodir, "/R/WaveQTL_preprocess_funcs.R"))

## Read phenotype data
phenoD = read.table(paste0(wd.path, "/phenoD.", ss), as.is = TRUE)
meanR.thresh = 2
res = WaveQTL_preprocess(Data = as.matrix(phenoD), library.read.depth =NULL , Covariates = NULL, meanR.thresh = meanR.thresh, no.QT = TRUE)
filteredWCs = res$filtered.WCs
norm.DNase = res$WCs

##save normaized data and useWCs information in output.path
this.path = paste0(wd.path, "/DNase.", ss, ".txt")
write.table(norm.DNase, file = this.path, quote= FALSE, row.names = FALSE, col.names = FALSE)
this.path = paste0(wd.path, "/use.", ss, ".txt")
cat(filteredWCs, file = this.path)


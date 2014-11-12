#!/usr/bin/env Rscript

## Aim : This file contains Rscript to compute area under roc curves from two methods on simulated data with different sample sizes and library read depth
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+


## for 70 30
## see.logLR.R
## case.name = c("fullread.70ind.over", "fullread.30ind.over", "fullread.10ind.over")
## case.name = c("halfread.70ind.over", "halfread.30ind.over", "2fullread.10ind.over")

## for 10
## see.logLR.RD.10ind.R
## case.name = c("fullread.10ind.over", "2fullread.10ind.over", "4fullread.10ind.over")

## for 6
## see.logLR.RD.6ind.R
## case.name = c("fullread.6ind.over", "2fullread.6ind.over", "4fullread.6ind.over")

## for 4
## see.logLR.RD.4ind.R
## case.name = c("fullread.4ind.over", "2fullread.4ind.over", "4fullread.4ind.over")




###############################
## Save data as a matrix form
## 5 by 4
## 5 (70, 30, 10, 6, 4)
## 4 (a half, full, 2full, 4full
###############################

AuROC.wave = matrix(data = NA, nc = 4, nr = 5)
AuROC.ms = matrix(data = NA, nc = 4, nr = 5)

status = c(rep(0, 500), rep(1, 500))

#################
## 70
#################
ind.ix = 1

case.name = c("fullread.70ind.over", "halfread.70ind.over")


##### ms ######
ms.null = vector("list", length(case.name))
ms.alt = vector("list", length(case.name))


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[1], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[1], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[1]] = logLR.null
ms.alt[[1]] = logLR.alt



load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[2], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[2], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res





sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[2]] = logLR.null
ms.alt[[2]] = logLR.alt


##### wave ######


wave.null = vector("list", length(case.name))
wave.alt = vector("list", length(case.name))


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[1], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[1], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[1]] = logLR.null
wave.alt[[1]] = logLR.alt


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[2], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[2], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res




sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[2]] = logLR.null
wave.alt[[2]] = logLR.alt


######################
## Compute AUC
######################

library("pROC")

ix = 1
logLR = c(ms.null[[ix]], ms.alt[[ix]])
AuROC.ms[ind.ix, 2] = auc(status, logLR)[1]

logLR = c(wave.null[[ix]], wave.alt[[ix]])
AuROC.wave[ind.ix, 2] = auc(status, logLR)[1]


ix = 2
logLR = c(ms.null[[ix]], ms.alt[[ix]])
AuROC.ms[ind.ix, 1] = auc(status, logLR)[1]

logLR = c(wave.null[[ix]], wave.alt[[ix]])
AuROC.wave[ind.ix, 1] = auc(status, logLR)[1]



#################
## 30
#################
ind.ix = 2

case.name = c("fullread.30ind.over", "halfread.30ind.over")


##### ms ######
ms.null = vector("list", length(case.name))
ms.alt = vector("list", length(case.name))


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[1], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[1], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[1]] = logLR.null
ms.alt[[1]] = logLR.alt



load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[2], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[2], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res





sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[2]] = logLR.null
ms.alt[[2]] = logLR.alt


##### wave ######


wave.null = vector("list", length(case.name))
wave.alt = vector("list", length(case.name))


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[1], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[1], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[1]] = logLR.null
wave.alt[[1]] = logLR.alt


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[2], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[2], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res




sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[2]] = logLR.null
wave.alt[[2]] = logLR.alt


######################
## Compute AUC
######################

#library("pROC")

ix = 1
logLR = c(ms.null[[ix]], ms.alt[[ix]])
AuROC.ms[ind.ix, 2] = auc(status, logLR)[1]

logLR = c(wave.null[[ix]], wave.alt[[ix]])
AuROC.wave[ind.ix, 2] = auc(status, logLR)[1]


ix = 2
logLR = c(ms.null[[ix]], ms.alt[[ix]])
AuROC.ms[ind.ix, 1] = auc(status, logLR)[1]

logLR = c(wave.null[[ix]], wave.alt[[ix]])
AuROC.wave[ind.ix, 1] = auc(status, logLR)[1]





#################
## 10
#################
ind.ix = 3


case.name = c("fullread.10ind.over", "2fullread.10ind.over", "4fullread.10ind.over")


##### ms ######
ms.null = vector("list", length(case.name))
ms.alt = vector("list", length(case.name))


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[1], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[1], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[1]] = logLR.null
ms.alt[[1]] = logLR.alt



load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[2], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[2], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res





sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[2]] = logLR.null
ms.alt[[2]] = logLR.alt



load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[3], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[3], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res





sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[3]] = logLR.null
ms.alt[[3]] = logLR.alt


##### wave ######


wave.null = vector("list", length(case.name))
wave.alt = vector("list", length(case.name))


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[1], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[1], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[1]] = logLR.null
wave.alt[[1]] = logLR.alt


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[2], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[2], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res




sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[2]] = logLR.null
wave.alt[[2]] = logLR.alt



load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[3], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[3], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res




sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[3]] = logLR.null
wave.alt[[3]] = logLR.alt




######################
## Compute AUC
######################

#library("pROC")

ix = 1
logLR = c(ms.null[[ix]], ms.alt[[ix]])
AuROC.ms[ind.ix, 2] = auc(status, logLR)[1]

logLR = c(wave.null[[ix]], wave.alt[[ix]])
AuROC.wave[ind.ix, 2] = auc(status, logLR)[1]


ix = 2
logLR = c(ms.null[[ix]], ms.alt[[ix]])
AuROC.ms[ind.ix, 3] = auc(status, logLR)[1]

logLR = c(wave.null[[ix]], wave.alt[[ix]])
AuROC.wave[ind.ix, 3] = auc(status, logLR)[1]


ix = 3
logLR = c(ms.null[[ix]], ms.alt[[ix]])
AuROC.ms[ind.ix, 4] = auc(status, logLR)[1]

logLR = c(wave.null[[ix]], wave.alt[[ix]])
AuROC.wave[ind.ix, 4] = auc(status, logLR)[1]




#################
## 6
#################
ind.ix = 4


case.name = c("fullread.6ind.over", "2fullread.6ind.over", "4fullread.6ind.over")


##### ms ######
ms.null = vector("list", length(case.name))
ms.alt = vector("list", length(case.name))


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[1], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[1], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[1]] = logLR.null
ms.alt[[1]] = logLR.alt



load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[2], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[2], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res





sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[2]] = logLR.null
ms.alt[[2]] = logLR.alt



load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[3], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[3], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res





sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[3]] = logLR.null
ms.alt[[3]] = logLR.alt


##### wave ######


wave.null = vector("list", length(case.name))
wave.alt = vector("list", length(case.name))


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[1], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[1], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[1]] = logLR.null
wave.alt[[1]] = logLR.alt


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[2], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[2], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res




sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[2]] = logLR.null
wave.alt[[2]] = logLR.alt



load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[3], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[3], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res




sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[3]] = logLR.null
wave.alt[[3]] = logLR.alt




######################
## Compute AUC
######################

#library("pROC")

ix = 1
logLR = c(ms.null[[ix]], ms.alt[[ix]])
AuROC.ms[ind.ix, 2] = auc(status, logLR)[1]

logLR = c(wave.null[[ix]], wave.alt[[ix]])
AuROC.wave[ind.ix, 2] = auc(status, logLR)[1]


ix = 2
logLR = c(ms.null[[ix]], ms.alt[[ix]])
AuROC.ms[ind.ix, 3] = auc(status, logLR)[1]

logLR = c(wave.null[[ix]], wave.alt[[ix]])
AuROC.wave[ind.ix, 3] = auc(status, logLR)[1]


ix = 3
logLR = c(ms.null[[ix]], ms.alt[[ix]])
AuROC.ms[ind.ix, 4] = auc(status, logLR)[1]

logLR = c(wave.null[[ix]], wave.alt[[ix]])
AuROC.wave[ind.ix, 4] = auc(status, logLR)[1]




#################
## 4
#################
ind.ix = 5


case.name = c("fullread.4ind.over", "2fullread.4ind.over", "4fullread.4ind.over")


##### ms ######
ms.null = vector("list", length(case.name))
ms.alt = vector("list", length(case.name))


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[1], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[1], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[1]] = logLR.null
ms.alt[[1]] = logLR.alt



load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[2], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[2], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res





sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[2]] = logLR.null
ms.alt[[2]] = logLR.alt



load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[3], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[3], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res





sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[3]] = logLR.null
ms.alt[[3]] = logLR.alt


##### wave ######


wave.null = vector("list", length(case.name))
wave.alt = vector("list", length(case.name))


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[1], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[1], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[1]] = logLR.null
wave.alt[[1]] = logLR.alt


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[2], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[2], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res




sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[2]] = logLR.null
wave.alt[[2]] = logLR.alt



load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[3], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[3], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res




sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[3]] = logLR.null
wave.alt[[3]] = logLR.alt




######################
## Compute AUC
######################

#library("pROC")

ix = 1
logLR = c(ms.null[[ix]], ms.alt[[ix]])
AuROC.ms[ind.ix, 2] = auc(status, logLR)[1]

logLR = c(wave.null[[ix]], wave.alt[[ix]])
AuROC.wave[ind.ix, 2] = auc(status, logLR)[1]


ix = 2
logLR = c(ms.null[[ix]], ms.alt[[ix]])
AuROC.ms[ind.ix, 3] = auc(status, logLR)[1]

logLR = c(wave.null[[ix]], wave.alt[[ix]])
AuROC.wave[ind.ix, 3] = auc(status, logLR)[1]


ix = 3
logLR = c(ms.null[[ix]], ms.alt[[ix]])
AuROC.ms[ind.ix, 4] = auc(status, logLR)[1]

logLR = c(wave.null[[ix]], wave.alt[[ix]])
AuROC.wave[ind.ix, 4] = auc(status, logLR)[1]



output.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/res/AuROC.wave"
write.table(AuROC.wave, file=output.path, quote=FALSE, row.names = c(70,30,10,6,4), col.names = FALSE)



output.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/res/AuROC.ms"
write.table(AuROC.ms, file=output.path, quote=FALSE, row.names = c(70,30,10,6,4), col.names = FALSE)




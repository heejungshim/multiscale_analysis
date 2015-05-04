## `preprocess.wave.run.multiscale.simulation.R' simulates data from Poisson distribution for wavelet analysis and multiseq, preprocess data for WaveQTL, and run multiseq on simulated data. 
##
##
## Example Usage : 
## /data/tools/R-3.0.3/bin/R CMD BATCH --no-save --no-restore "--args seed=$SGE_TASK_ID geno.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_simple_v1/data/geno70.dat' sig0.path='~/multiscale_analysis/analysis/simulation/sample_size/simulation_simple_v1/data/alt.sig0' sig1.path='~/multiscale_analysis/analysis/simulation/sample_size/simulation_simple_v1/data/alt.sig1' wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_simple_v1/alt/' output.dir.name='fullread.10ind' multipleSig=1" /mnt/lustre/home/shim/multiscale_analysis/src/R/preprocess.wave.run.multiscale.simulation.R
## 
##
##
## Copyright (C) 2014 Heejung Shim
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.




setwd("/mnt/lustre/home/shim/multiscale_analysis")

library("multiseq")
library("ashr")

multiscale.analysis.repodir <- scan(".multiscale_analysis.repodir.txt", what=character())
source(paste0(multiscale.analysis.repodir, "/src/R/my.utils.R"))

WaveQTL.repodir <- scan(".WaveQTL.repodir.txt", what=character())


##seed=1
##geno.path='~/multiscale_analysis/analysis/simulation/sample_size/simulation_simple_v1/data/geno10.dat'
##sig0.path='~/multiscale_analysis/analysis/simulation/sample_size/simulation_simple_v1/data/alt.sig0'
##sig1.path='~/multiscale_analysis/analysis/simulation/sample_size/simulation_simple_v1/data/alt.sig1'
##wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_simple_v1/alt/'
##output.dir.name='fullread.10ind'
##multipleSig=1 


args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))
eval(parse(text=args[[2]]))
eval(parse(text=args[[3]]))
eval(parse(text=args[[4]]))
eval(parse(text=args[[5]]))
eval(parse(text=args[[6]]))
if(length(args) > 6){
  eval(parse(text=args[[7]]))
  if(multipleSig == 1){
    sig0.path = paste0(sig0.path, ".", seed)
    sig1.path = paste0(sig1.path, ".", seed)
  }
}

setwd(wd.path)



#####################
# sample data
#####################

# read genotype data
genoD = round(as.numeric(scan(geno.path, what=double())))
numSam = length(genoD)

# read signal
mu0 = scan(file = sig0.path, what = double())
mu1 = scan(file = sig1.path, what = double())
numBPs = length(mu0)

# phenotype data
phenoD = matrix(data=NA, nr= length(genoD), nc = numBPs)

# let's sample!!!
set.seed(seed)

# geno = 0
wh0 = which(genoD == 0)
if(length(wh0) > 0){
  numSam = length(wh0)
  phenoD[wh0,] = matrix(data = rpois(numSam*numBPs, mu0), byrow = T, nr = numSam)
}
# geno = 1
wh1 = which(genoD == 1)
if(length(wh1) > 0){
  numSam = length(wh1)
  phenoD[wh1,] = matrix(data = rpois(numSam*numBPs, mu1), byrow = T, nr = numSam)
}


#########################################
# data preprocessing for wavelet analysis
#########################################

source(paste0(WaveQTL.repodir, "/R/WaveQTL_preprocess_funcs.R"))

meanR.thresh = 2
res = WaveQTL_preprocess(Data = phenoD, library.read.depth =NULL , Covariates = NULL, meanR.thresh = meanR.thresh)

filteredWCs = res$filtered.WCs
norm.DNase = res$WCs

##save normaized data and useWCs information in output.path
out.dir.path = paste0(wd.path, "wave/", output.dir.name, ".data/") 
if(!file.exists(out.dir.path)){
  dir.create(out.dir.path)
}

this.path = paste0(out.dir.path, "DNase.", seed, ".txt")
write.table(norm.DNase, file = this.path, quote= FALSE, row.names = FALSE, col.names = FALSE)
this.path = paste0(out.dir.path, "use.", seed, ".txt")
cat(filteredWCs, file = this.path)


numPerm = NULL
numSig = 10
# perform test 
res = permutation.logLR(pheno.dat = phenoD, geno.dat = genoD, library.read.depth = NULL, numPerm = numPerm, numSig= numSig, use.default.compute.logLR = TRUE, cxx=TRUE)



# write output
out.dir.path = paste0(wd.path, "multiscale/", output.dir.name, ".output") 
if(!file.exists(out.dir.path)){
    dir.create(out.dir.path)
}


write.table(res$logLR, file = paste0(out.dir.path, "/res.", seed, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE)





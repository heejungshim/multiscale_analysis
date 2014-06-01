## `run.multiscale.on.simulation.R' simulates data by thinning and perform
## permutation-based test using a multiscale Poisson model. See the input arguments. 
## Example Usage : R CMD BATCH --no-save --no-restore "--args seed=$SGE_TASK_ID numPerm=1000 numSig=10 geno.path='/mnt/lustre/home/shim/wavelets/revision/analysis/simulation_578/data/geno.$SGE_TASK_ID.dat' raw.dat.path='/mnt/lustre/home/shim/wavelets/revision/analysis/simulation_578/data/raw.$SGE_TASK_ID.dat' ratio.path='/mnt/lustre/home/shim/wavelets/revision/analysis/simulation_578/data/smooth.pro1.21.$SGE_TASK_ID' scale.level=0.8 wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/multiscale/' read.depth.ratio=NULL" /mnt/lustre/home/shim/multiscale_analysis/src/R/run.multiscale.on.simulation.R
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


#seed=1
#numPerm=20
#numSig=3
#geno.path='/mnt/lustre/home/shim/wavelets/revision/analysis/simulation_578/data/geno.1.dat'
#raw.dat.path='/mnt/lustre/home/shim/wavelets/revision/analysis/simulation_578/data/raw.1.dat'
#ratio.path='/mnt/lustre/home/shim/wavelets/revision/analysis/simulation_578/data/smooth.pro1.21.1'
#scale.level=NULL
#wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/multiscale/'
#read.depth.ratio=0.5

args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))
eval(parse(text=args[[2]]))
eval(parse(text=args[[3]]))
eval(parse(text=args[[4]]))
eval(parse(text=args[[5]]))
eval(parse(text=args[[6]]))
eval(parse(text=args[[7]]))
eval(parse(text=args[[8]]))
eval(parse(text=args[[9]]))


 

setwd(wd.path)






#####################
# sample data
#####################

# read raw data from which we will sample
raw.data = read.table(raw.dat.path, as.is = TRUE)
raw.data.T = as.numeric(apply(raw.data, 2, sum))

# change read depth
if(!is.null(read.depth.ratio)){
    raw.data.T = floor(raw.data.T*read.depth.ratio)
}

mu0.sig = rep(1/70, 1024)
if(is.null(ratio.path)){
    mu1.sig = mu0.sig
}else{
    ratio = as.numeric(scan(ratio.path, what=double()))
    if(!is.null(scale.level)){
        ratio = 1 + scale.level*(ratio - 1)
    }
    mu1.sig = mu0.sig*ratio
}
    

# ok!!
# upper and lower bound!
trunc.fun = function(x){
    x = max(0, x)
    return(min(1,x))
}
mu0.sig = sapply(mu0.sig, trunc.fun)
mu1.sig = sapply(mu1.sig, trunc.fun)

# read genotype data
genoD = round(as.numeric(scan(geno.path, what=double())))

# phenotype data
phenoD = matrix(data=NA, nr= length(genoD), nc = length(raw.data.T))

# let's sample!!!
set.seed(seed)
# geno = 0
wh0 = which(genoD == 0)
if(length(wh0) > 0){
    dat.T = rbinom(length(raw.data.T)*length(wh0), raw.data.T, mu0.sig) 
    phenoD[wh0,] = matrix(data=dat.T, nr = length(wh0), byrow = TRUE)
}  
# geno = 1
wh1 = which(genoD == 1)
if(length(wh1) > 0){
    dat.T = rbinom(length(raw.data.T)*length(wh1), raw.data.T, mu1.sig) 
    phenoD[wh1,] = matrix(data=dat.T, nr = length(wh1), byrow = TRUE)
}  



# perform test 
res2 = permutation.logLR(pheno.dat = phenoD, geno.dat = genoD, library.read.depth = NULL, numPerm = numPerm, numSig= numSig, use.default.compute.logLR = TRUE, cxx=TRUE)


# write output
out.dir.path = paste0(wd.path, "output") 
if(!file.exists(out.dir.path)){
    dir.create(out.dir.path)
}

write.table(res$logLR, file = paste0(out.dir.path, "/res.", seed, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE)

pval.path = paste0(out.dir.path, "/pval.", seed, ".out")
cat(1, file =pval.path)
cat("\n", file = pval.path, append = TRUE)
cat(res$Count_stop, file = pval.path, append = TRUE)
cat("\n", file =pval.path, append = TRUE)
cat(res$Count_sig, file = pval.path, append = TRUE)
cat("\n", file = pval.path, append = TRUE)	
cat(res$pval, file = pval.path, append = TRUE)



out.dir.path = paste0(wd.path, "warnings") 
if(!file.exists(out.dir.path)){
    dir.create(out.dir.path)
}


warn.path = paste0(out.dir.path, "/warnings.", seed, ".out")
cat(warnings(), file =warn.path)

















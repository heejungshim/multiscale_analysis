## `preprocess.wave.window.simulation.R' simulates data by thinning and pre-process data for wavelet analysis (if wavelet.preprocess== TRUE) and/or for window based analysis (if window.preprocess==TRUE).
##
##
## Example Usage : 
## /data/tools/R-3.0.3/bin/R CMD BATCH --no-save --no-restore "--args seed=$SGE_TASK_ID geno.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/data/geno70.dat' raw.dat.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/data/raw.dat' ratio.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/data/smooth.ratio' scale.level=0.5 wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/' read.depth.ratio=NULL output.dir.name='fullread.70ind' wavelet.preprocess=TRUE window.preprocess=TRUE over.dispersion=NULL" /mnt/lustre/home/shim/multiscale_analysis/src/R/preprocess.wave.window.simulation.R
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




#seed=1
#geno.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/data/geno10.dat'
#raw.dat.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/data/raw.dat'
#ratio.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/data/smooth.ratio'
#scale.level=0.5
#wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/'
#read.depth.ratio=NULL
#output.dir.name='fullread.10ind'
#wavelet.preprocess=TRUE
#window.preprocess=TRUE
#over.dispersion=1/70/70/10


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
eval(parse(text=args[[10]]))
eval(parse(text=args[[11]]))


 

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
    phenoD[wh0,] = sample.from.Binomial.with.Overdispersion(num.sam = length(wh0), total.count = raw.data.T, mu.sig = mu0.sig, over.dispersion = over.dispersion)
}  
# geno = 1
wh1 = which(genoD == 1)
if(length(wh1) > 0){
    phenoD[wh1,] = sample.from.Binomial.with.Overdispersion(num.sam = length(wh1), total.count = raw.data.T, mu.sig = mu1.sig, over.dispersion = over.dispersion)
}




#########################################
# data preprocessing for wavelet analysis
#########################################

if(wavelet.preprocess){
    source("/mnt/lustre/home/shim/BIMBAM_wavelets/WaveQTL_example/Wave_preprocess.R")
    meanR.thresh = 2
    res = Wave_preprocess(Data = phenoD, Read.depth =NULL , C = NULL, meanR.thresh = meanR.thresh)

    filteredWCs = res$filtered.WCs
    norm.DNase = res$WCs

    #save normaized data and useWCs information in output.path
    out.dir.path = paste0(wd.path, "wave/", output.dir.name, ".data/") 
    if(!file.exists(out.dir.path)){
        dir.create(out.dir.path)
    }

    this.path = paste0(out.dir.path, "DNase.", seed, ".txt")
    write.table(norm.DNase, file = this.path, quote= FALSE, row.names = FALSE, col.names = FALSE)
    this.path = paste0(out.dir.path, "use.", seed, ".txt")
    cat(filteredWCs, file = this.path)

}


#########################################
# data preprocessing for window based approach 
# modify a code from
# /mnt/lustre/home/shim/wavelets/data/DNase/command_01_step1/for_sel_data_new_Jack_01_step1/ready_100window.R
# /mnt/lustre/home/shim/wavelets/data/DNase/command_01_step1/for_sel_data_new_Jack_normal_01_step1/QN2.R
#########################################

if(window.preprocess){

    numBPs = dim(phenoD)[2]
    numC = numBPs%/%100
    numIND = dim(phenoD)[1]

    mat = matrix(data=NA, nc = numC, nr = numIND)
    st = 1
    for(c in 1:(numC-1)){
        en = st + 100 - 1
        den = en - st + 1
        if(den > 0){
            mat[,c] = apply(phenoD[,st:en], 1, sum)/den
        }else{
            mat[,c] = 0
        }
        st = en + 1
    }
    en = numBPs

    den = en - st + 1
    if(den > 0){
        mat[,numC] = apply(phenoD[,st:en], 1, sum)/den
    }else{
        mat[,numC] = 0
    }

    QT_randomTie <- function(x) {
        x.rank = rank(x, ties.method="random")
        return(qqnorm(x.rank,plot.it = F)$x)
    }

    QT_dat = apply(mat, 2, QT_randomTie)

    # save QT data
    out.dir.path = paste0(wd.path, "window/", output.dir.name, ".data/") 
    if(!file.exists(out.dir.path)){
        dir.create(out.dir.path)
    }

    this.path = paste0(out.dir.path, "JAR.", seed, ".txt")
    write.table(QT_dat, file = this.path, quote= FALSE, row.names = FALSE, col.names = FALSE)
}



























setwd("/mnt/lustre/home/shim/multiscale_analysis")

library("multiseq")
library("ashr")


multiscale.analysis.repodir <- scan(".multiscale_analysis.repodir.txt", what=character())
source(paste0(multiscale.analysis.repodir, "/src/R/my.utils.R"))



seed=89
numPerm=1000
numSig=10
geno.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/data/geno70.dat'
raw.dat.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/data/raw.dat.89'
ratio.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/data/smooth.ratio.89'
scale.level=NULL
wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/multiscale/'
read.depth.ratio=NULL
output.dir.name='fullread.70ind'
over.dispersion=NULL




setwd(wd.path)



### create file for warning message
#out.dir.path = paste0(wd.path, output.dir.name, ".warnings") 
#if(!file.exists(out.dir.path)){
#    dir.create(out.dir.path)
#}

#warn.path = paste0(out.dir.path, "/warnings.", seed, ".txt")
#warnings.file <- file(warn.path, open="wt")
#sink(warnings.file, type="message")



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
#Warning message:
#In rbinom(length(mu.sig) * num.sam, total.count, mu.sig) : NAs produced


# geno = 1
wh1 = which(genoD == 1)
if(length(wh1) > 0){
    phenoD[wh1,] = sample.from.Binomial.with.Overdispersion(num.sam = length(wh1), total.count = raw.data.T, mu.sig = mu1.sig, over.dispersion = over.dispersion)
}
#Warning message:
#In rbinom(length(mu.sig) * num.sam, total.count, mu.sig) : NAs produced


## let's look at what's going on inside of function


    phenoD[wh1,] = sample.from.Binomial.with.Overdispersion(num.sam = length(wh1), total.count = raw.data.T, mu.sig = mu1.sig, over.dispersion = over.dispersion)



num.sam = length(wh1)
total.count = raw.data.T
mu.sig = mu1.sig
over.dispersion = over.dispersion


rbinom(1, total.count[829], mu.sig[829])
#[1] NaN
#Warning message:
#In rbinom(1, total.count[829], mu.sig[829]) : NAs produced
total.count[829]
#[1] 0.5
mu.sig[829]
#[1] 0.01828146


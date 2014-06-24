


setwd("/mnt/lustre/home/shim/multiscale_analysis")

library("multiseq")
library("ashr")

multiscale.analysis.repodir <- scan(".multiscale_analysis.repodir.txt", what=character())
#source(paste0(multiscale.analysis.repodir, "/src/R/my.utils.R"))
multiseq.repodir <- scan(".multiseq.repodir.txt", what=character())
source(paste0(multiseq.repodir, "/package/multiseq/R/multiseq.R"))

ash.repodir <- scan(".ash.repodir.txt", what=character())
source(paste0(ash.repodir, "/package/ashr/R/ash.R"))
source(paste0(ash.repodir, "/package/ashr/R/mix.R"))



setwd("/mnt/lustre/home/shim/multiscale_analysis/debug/bin_vs_multi/code/")

load("/mnt/lustre/home/shim/multiscale_analysis/debug/bin_vs_multi/data/sample_data_binmulti_large.Robj")
 
locus.end
#[1] 25821184
locus.start
#[1] 25690113
peak.start
#[1] 25799098
peak.end
#[1] 25799441
read.depth
#[1] 16335812 18197248 24225586 12378544
2^17
#[1] 131072
locus.end - locus.start + 1
#[1] 131072
dim(sim.data.bin.alt[[1]])
#[1]      4 131072
dim(sim.data.bin.null[[1]])
#[1]      4 131072
dim(sim.data.multi.alt[[1]])
#[1]      4 131072
dim(sim.data.multi.null[[1]])
#[1]      4 131072
dim(sim.data.nb.alt[[1]])
#[1]      4 131072
dim(sim.data.nb.null[[1]])
#[1]      4 131072

bin.alt = multiseq(x = sim.data.bin.alt[[1]], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

bin.null = multiseq(x = sim.data.bin.null[[1]], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

multi.alt = multiseq(x = sim.data.multi.alt[[1]], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)
#Warning message:
#In EMest(betahat[completeobs], lambda1 * sebetahat[completeobs] +  :
#  EM algorithm in function mixEM failed to converge. Results may be unreliable. Try increasing maxiter and rerunning.

multi.null = multiseq(x = sim.data.multi.null[[1]], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)


bin.alt$logLR$value
#[1] 0.01503283
bin.null$logLR$value
#[1] 0.02489572
multi.alt$logLR$value
#[1] 7.351711
multi.null$logLR$value
#[1] 0.1439083


round(bin.null$logLR$scales,3)
# [1] 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
#[13] 0.000 0.000 0.000 0.000 0.000 0.025
round(bin.alt$logLR$scales,3)
# [1] 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.015 0.000
#[13] 0.000 0.000 0.000 0.000 0.000 0.000
round(multi.null$logLR$scales,3)
# [1] 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
#[13] 0.000 0.000 0.089 0.000 0.054 0.000
round(multi.alt$logLR$scales,3)
# [1] 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.318 2.038 0.719
#[13] 1.889 2.371 0.016 0.000 0.000 0.000


round(multi.alt$fitted.g[[14]]$pi,2)
#[1] 0.29 0.00 0.00 0.00 0.00 0.39 0.32 0.00 0.00
round(multi.alt$fitted.g[[14]]$sd,2)
#[1] 0.00 0.01 0.03 0.06 0.11 0.22 0.45 0.90 1.80
round(bin.alt$fitted.g[[14]]$pi,2)
#[1] 1 0 0 0 0 0 0 0 0
round(bin.alt$fitted.g[[14]]$sd,2)
#[1] 0.00 0.01 0.02 0.04 0.08 0.16 0.32 0.64 1.28


round(multi.null$fitted.g[[14]]$pi,2)
# [1] 1 0 0 0 0 0 0 0 0 0
round(multi.null$fitted.g[[14]]$sd,2)
# [1] 0.00 0.01 0.02 0.03 0.06 0.12 0.25 0.50 1.00 1.99
round(bin.null$fitted.g[[14]]$pi,2)
#[1] 1 0 0 0 0 0 0 0 0
round(bin.null$fitted.g[[14]]$sd,2)
#[1] 0.00 0.01 0.02 0.03 0.07 0.14 0.27 0.54 1.08


round(multi.alt$fitted.g[[13]]$pi,2)
# [1] 0.01 0.00 0.00 0.00 0.64 0.02 0.30 0.03 0.00 0.00
round(multi.alt$fitted.g[[13]]$sd,2)
# [1] 0.00 0.01 0.02 0.05 0.10 0.20 0.39 0.78 1.56 3.12
round(bin.alt$fitted.g[[13]]$pi,2)
# [1] 1 0 0 0 0 0 0 0 0 0
round(bin.alt$fitted.g[[13]]$sd,2)
# [1] 0.00 0.01 0.02 0.04 0.09 0.18 0.35 0.71 1.42 2.83
 
round(multi.null$fitted.g[[13]]$pi,2)
# [1] 1 0 0 0 0 0 0 0 0 0
round(multi.null$fitted.g[[13]]$sd,2)
# [1] 0.00 0.01 0.02 0.04 0.09 0.18 0.35 0.71 1.42 2.84
round(bin.null$fitted.g[[13]]$pi,2)
# [1] 1 0 0 0 0 0 0 0 0 0
round(bin.null$fitted.g[[13]]$sd,2)
# [1] 0.00 0.01 0.03 0.05 0.10 0.21 0.42 0.83 1.66 3.32


round(multi.alt$fitted.g[[11]]$pi,2)
# [1] 0.98 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.02 0.00 0.00
round(multi.alt$fitted.g[[11]]$sd,2)
# [1] 0.00 0.02 0.03 0.06 0.12 0.24 0.49 0.98 1.95 3.91 7.81
round(bin.alt$fitted.g[[11]]$pi,2)
# [1] 0.99 0.00 0.00 0.00 0.00 0.00 0.00 0.01 0.00 0.00
round(bin.alt$fitted.g[[11]]$sd,2)
# [1] 0.00 0.02 0.03 0.06 0.12 0.25 0.50 0.99 1.99 3.97


round(multi.null$fitted.g[[11]]$pi,2)
# [1] 1 0 0 0 0 0 0 0 0 0 0
round(multi.null$fitted.g[[11]]$sd,2)
# [1] 0.00 0.01 0.03 0.05 0.10 0.21 0.41 0.83 1.66 3.31 6.63
round(bin.null$fitted.g[[11]]$pi,2)
# [1] 1 0 0 0 0 0 0 0 0 0 0
round(bin.null$fitted.g[[11]]$sd,2)
# [1] 0.00 0.02 0.03 0.06 0.12 0.25 0.49 0.98 1.97 3.94 7.87


round(multi.alt$fitted.g[[17]]$pi,2)
#[1] 1 0 0 0 0 0 0 0
round(multi.alt$fitted.g[[17]]$sd,2)
#[1] 0.00 0.01 0.01 0.02 0.05 0.09 0.19 0.38
round(bin.alt$fitted.g[[17]]$pi,2)
#[1] 1 0 0 0 0 0 0
round(bin.alt$fitted.g[[17]]$sd,2)
#[1] 0.00 0.01 0.01 0.03 0.05 0.10 0.20


round(multi.null$fitted.g[[17]]$pi,2)
#[1] 0 0 0 0 1 0 0 0
round(multi.null$fitted.g[[17]]$sd,2)
#[1] 0.00 0.01 0.01 0.02 0.05 0.09 0.18 0.36
round(bin.null$fitted.g[[17]]$pi,2)
#[1] 1 0 0 0 0 0 0
round(bin.null$fitted.g[[17]]$sd,2)
#[1] 0.00 0.00 0.01 0.02 0.04 0.08 0.15


round(multi.alt$fitted.g[[18]]$pi,2)
#[1] 1 0 0 0 0
round(multi.alt$fitted.g[[18]]$sd,2)
#[1] 0.00 0.00 0.01 0.01 0.03
round(bin.alt$fitted.g[[18]]$pi,2)
#[1] 1 0 0 0 0
round(bin.alt$fitted.g[[18]]$sd,2)
#[1] 0.00 0.01 0.01 0.02 0.04


round(multi.null$fitted.g[[18]]$pi,2)
#[1] 1 0 0 0 0
round(multi.null$fitted.g[[18]]$sd,2)
#[1] 0.00 0.00 0.01 0.01 0.03
round(bin.null$fitted.g[[18]]$pi,2)
#[1] 0 0 0 0 1 0
round(bin.null$fitted.g[[18]]$sd,2)
#[1] 0.00 0.00 0.01 0.01 0.02 0.04





######## Look at what's going on inside of multiseq


#multiseq = function(x,g=NULL,read.depth = NULL,reflect=FALSE,baseline="inter",minobs=1,pseudocounts=0.5,all=FALSE,center=FALSE,repara=TRUE,forcebin=FALSE,lm.approx=TRUE,disp=c("add","mult"),nullcheck=TRUE,pointmass=TRUE,prior="nullbiased",gridmult=2,mixsd=NULL,VB=FALSE,shape.eff=FALSE,cxx=TRUE, computelogLR = FALSE, maxlogLR = NULL){


#multi.null = multiseq( x = sim.data.multi.null[[1]], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

g = c(0, 0, 1, 1)
read.depth = read.depth
prior = "uniform"
onlylogLR = TRUE

reflect=FALSE
baseline="inter"
minobs=1
pseudocounts=0.5
all=FALSE
center=FALSE
repara=TRUE
forcebin=FALSE
lm.approx=TRUE
disp="add"
nullcheck=TRUE
pointmass=TRUE
gridmult=2
mixsd=NULL
VB=FALSE
shape.eff=FALSE
cxx=TRUE
smoothing=TRUE
cyclespin=TRUE
reverse=TRUE
maxlogLR=NULL
set.fitted.g=NULL
set.fitted.g.intercept=NULL
get.fitted.g=TRUE
listy=NULL






    if(onlylogLR){
        if(is.null(g)) stop("Error: g should be provided to compute logLR (onlylogLR = TRUE)")
        if(pointmass != TRUE) stop("Error: logLR can be computed only when pointmass = TRUE")
        reverse=FALSE
    }
    if (!smoothing) reverse = FALSE
    if (!cyclespin) {reverse = FALSE; warning("Reversing wavelet not implemented here when cyclespin=FALSE, setting reverse=FALSE")}
    #to do: check other input parameters



nsig = 4
J = 17
n = 2^17

if(!is.null(g)){
        if(is.factor(g))
            g.num = as.numeric(levels(g))[g]
        else{
            g.num = g
            if(length(unique(g)) == 2)
                g = factor(g)
        }
    }


## run for multinomial alt data

x = sim.data.multi.alt[[1]]


    fitted.g=list()
    fitted.g.intercept=list()
    logLR=NULL
    sumlogLR=NULL
    finite.logLR=NULL





xRowSums = rowSums(x)
y.rate = matrix(c(xRowSums, read.depth-xRowSums), ncol=2)
zdat.rate = as.vector(glm.approx(y.rate, g=g, center=center, repara=repara, lm.approx=lm.approx, disp=disp))
zdat.rate.ash = ash(zdat.rate[3], zdat.rate[4], prior=prior, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=set.fitted.g[[J+1]])



fitted.g[[J+1]] = zdat.rate.ash$fitted.g
logLR[J+1] = zdat.rate.ash$logLR

y = matrix(nrow=nsig, ncol=2*J*n); for(i in 1:nsig){tt = ParentTItable(x[i,]); y[i,] = as.vector(t(tt$parent))}

zdat = glm.approx(y, g, minobs=minobs, pseudocounts=pseudocounts, center=center, all=all, forcebin=forcebin, repara=repara, lm.approx=lm.approx, disp=disp)


for(j in 1:J){

        spins = 2^j
        ind = ((j-1)*n+1):(j*n)
        zdat.ash = ash(zdat[3,ind], zdat[4,ind], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=set.fitted.g[[j]])
        
        fitted.g[[j]] = zdat.ash$fitted.g
        
        logLR[j] = zdat.ash$logLR/spins
}


sumlogLR = sum(logLR) # combine logLR from different scales



y.multi.alt = y
zdat.multi.alt = zdat
fitted.g.multi.alt =  fitted.g
logLR.multi.alt = logLR


# run on binomial alternative data

x = sim.data.bin.alt[[1]]

    fitted.g=list()
    fitted.g.intercept=list()
    logLR=NULL
    sumlogLR=NULL
    finite.logLR=NULL





xRowSums = rowSums(x)
y.rate = matrix(c(xRowSums, read.depth-xRowSums), ncol=2)
zdat.rate = as.vector(glm.approx(y.rate, g=g, center=center, repara=repara, lm.approx=lm.approx, disp=disp))
zdat.rate.ash = ash(zdat.rate[3], zdat.rate[4], prior=prior, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=set.fitted.g[[J+1]])



fitted.g[[J+1]] = zdat.rate.ash$fitted.g
logLR[J+1] = zdat.rate.ash$logLR

y = matrix(nrow=nsig, ncol=2*J*n); for(i in 1:nsig){tt = ParentTItable(x[i,]); y[i,] = as.vector(t(tt$parent))}

zdat = glm.approx(y, g, minobs=minobs, pseudocounts=pseudocounts, center=center, all=all, forcebin=forcebin, repara=repara, lm.approx=lm.approx, disp=disp)


for(j in 1:J){

        spins = 2^j
        ind = ((j-1)*n+1):(j*n)
        zdat.ash = ash(zdat[3,ind], zdat[4,ind], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=set.fitted.g[[j]])
        
        fitted.g[[j]] = zdat.ash$fitted.g
        
        logLR[j] = zdat.ash$logLR/spins
}


sumlogLR = sum(logLR) # combine logLR from different scales




y.bin.alt = y
zdat.bin.alt = zdat
fitted.g.bin.alt =  fitted.g
logLR.bin.alt = logLR


round(logLR.bin.alt,2)
# [1] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.02 0.00 0.00 0.00 0.00
#[16] 0.00 0.00 0.00
round(logLR.multi.alt,2)
# [1] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.32 2.04 0.72 1.89 2.37 0.02
#[16] 0.00 0.00 0.00
logLR.multi.alt[c(11,13,14)]
#[1] 2.037905 1.889064 2.371159
logLR.bin.alt[c(11,13,14)]
#[1] 1.502967e-02 3.162792e-06 0.000000e+00



### let's check if there is NA issue! I don't think so..

j = 11
ind = ((j-1)*n+1):(j*n)
sum(is.na(zdat.multi.alt[4,ind]))
# 3673
sum(is.na(zdat.bin.alt[4,ind]))
# 3742
length(ind)
# 131072

j = 13
ind = ((j-1)*n+1):(j*n)
sum(is.na(zdat.multi.alt[3,ind]))
# 0
sum(is.na(zdat.bin.alt[3,ind]))
# 0 
length(ind)

j = 14
ind = ((j-1)*n+1):(j*n)
sum(is.na(zdat.multi.alt[3,ind]))
# 0 
sum(is.na(zdat.bin.alt[3,ind]))
# 0 
length(ind)
# 131072


### let's check what happens if we use the other case's pi and sd
### without cyclepsin???


j = 14
ind = ((j-1)*n+1):(j*n)
length(ind)/(2^j)
st = min(ind)
en = st + 8 - 1

# binomial

zdat.ash.bin.alt = ash(zdat.bin.alt[3,st:en], zdat.bin.alt[4,st:en], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=set.fitted.g[[j]])

zdat.ash.bin.alt$logLR
#[1] 0.0006606527

zdat.ash.bin.alt = ash(zdat.bin.alt[3,ind], zdat.bin.alt[4,ind], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=set.fitted.g[[j]])

zdat.ash.bin.alt$logLR/(2^j)
#[1] 0


# multinomail

zdat.ash.multi.alt = ash(zdat.multi.alt[3,st:en], zdat.multi.alt[4,st:en], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=set.fitted.g[[j]])

zdat.ash.multi.alt$logLR
# [1] 6.611666

zdat.ash.multi.alt = ash(zdat.multi.alt[3,ind], zdat.multi.alt[4,ind], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR)

zdat.ash.multi.alt$logLR/(2^j)
# [1] 2.371159


# use the other data set's pi and sd
# bin <- use multi's pi and sd


zdat.ash.bin.alt = ash(zdat.bin.alt[3,st:en], zdat.bin.alt[4,st:en], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=fitted.g.multi.alt[[j]])

zdat.ash.bin.alt$logLR
#[1] 0

zdat.ash.bin.alt = ash(zdat.bin.alt[3,ind], zdat.bin.alt[4,ind], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=fitted.g.multi.alt[[j]])

zdat.ash.bin.alt$logLR/(2^j)
#[1] 0

# multi <- use binomial's pi and sd

zdat.ash.multi.alt = ash(zdat.multi.alt[3,st:en], zdat.multi.alt[4,st:en], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=fitted.g.bin.alt[[j]])

zdat.ash.multi.alt$logLR
# [1] 0

zdat.ash.multi.alt = ash(zdat.multi.alt[3,ind], zdat.multi.alt[4,ind], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=fitted.g.bin.alt[[j]])

zdat.ash.multi.alt$logLR/(2^j)
# [1] 0




######### let's look at MLE and standard error of MLE

zdat.bin.alt[3, st:en]
#[1] -0.380649174  0.003458425  0.231750418  0.188903955  0.109014227
#[6] -0.082618808 -0.336995397 -0.078746507
zdat.bin.alt[4, st:en]
#[1] 0.3420968 0.3504416 0.3106235 0.2464266 0.1479994 0.2800948 0.1798911
#[8] 0.1873700
 
zdat.multi.alt[3, st:en]
#[1]  0.152572863 -0.002395917  0.183514703  0.404047936 -0.051065806
#[6] -0.101343029 -0.810564823 -0.424790109
zdat.multi.alt[4, st:en]
#[1] 0.6141899 0.2700118 0.3115320 0.2507136 0.2053155 0.2064659 0.1814805
#[8] 0.1886367


zdat.bin.alt[3,st:en]/zdat.bin.alt[4, st:en]
#[1] -1.112694333  0.009868762  0.746081311  0.766573033  0.736585495
#[6] -0.294967341 -1.873330354 -0.420272674
zdat.multi.alt[3,st:en]/zdat.multi.alt[4, st:en]
#[1]  0.248413167 -0.008873381  0.589071690  1.611591568 -0.248718681
#[6] -0.490846347 -4.466402488 -2.251895683



# without the last two observations

zdat.ash.bin.alt = ash(zdat.bin.alt[3,st:(en-2)], zdat.bin.alt[4,st:(en-2)], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=set.fitted.g[[j]])

zdat.ash.bin.alt$logLR
# 0

zdat.ash.multi.alt = ash(zdat.multi.alt[3,st:(en-2)], zdat.multi.alt[4,st:(en-2)], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=set.fitted.g[[j]])

zdat.ash.multi.alt$logLR
# 0


# without one big observation

zdat.ash.bin.alt = ash(zdat.bin.alt[3,(st:en)[-7]], zdat.bin.alt[4,(st:en)[-7]], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=set.fitted.g[[j]])

zdat.ash.bin.alt$logLR
# 0

zdat.ash.multi.alt = ash(zdat.multi.alt[3,(st:en)[-7]], zdat.multi.alt[4,(st:en)[-7]], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=set.fitted.g[[j]])

zdat.ash.multi.alt$logLR
# 0.3018...



## Use only the last two observations 

zdat.ash.bin.alt = ash(zdat.bin.alt[3,(st:en)[7:8]], zdat.bin.alt[4,(st:en)[7:8]], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=set.fitted.g[[j]])

zdat.ash.bin.alt$logLR
# 0.2371

zdat.ash.multi.alt = ash(zdat.multi.alt[3,(st:en)[7:8]], zdat.multi.alt[4,(st:en)[7:8]], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR, g=set.fitted.g[[j]])

zdat.ash.multi.alt$logLR
# 8.916604


###### let's look at real count data correpsonding to those MLE and standard of MLE

dim(y.bin.alt)
dim(zdat.bin.alt)
s.IX = ((1:(dim(zdat.bin.alt)[2]))*2)[st:en]
f.IX = ((1:(dim(zdat.bin.alt)[2]))*2)[st:en] - 1

y.bin.alt[,s.IX]
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#[1,]   12   33   16   41   63   34   54   52
#[2,]   17   31   19   55   77   47   56   43
#[3,]   21   38   27   72   87   58   91   69
#[4,]   19   25   10   30   44   26   58   38
y.bin.alt[,f.IX]
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#[1,]   14   23   20   19  123   74   56   60
#[2,]   20   34   21   24  147   65   65   67
#[3,]   20   42   33   35  177   80   73   83
#[4,]   12   16   22   20  105   50   44   49

y.multi.alt[,s.IX]
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#[1,]   24   24   11   46   53   35   52   43
#[2,]   15   28   25   62   76   33   53   45
#[3,]   20   43   25   62   93   59  117   85
#[4,]    8   23    9   40   46   27   57   42
y.multi.alt[,f.IX]
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#[1,]   15   18   18   19  138   56   65   62
#[2,]   26   27   27   19  136   61   67   62
#[3,]   18   38   35   35  202   80   59   85
#[4,]   14   19   16   19   84   54   38   32

round(y.bin.alt[,s.IX]/(y.bin.alt[,s.IX] + y.bin.alt[,f.IX]),2)
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#[1,] 0.46 0.59 0.44 0.68 0.34 0.31 0.49 0.46
#[2,] 0.46 0.48 0.48 0.70 0.34 0.42 0.46 0.39
#[3,] 0.51 0.48 0.45 0.67 0.33 0.42 0.55 0.45
#[4,] 0.61 0.61 0.31 0.60 0.30 0.34 0.57 0.44

round(y.multi.alt[,s.IX]/(y.multi.alt[,s.IX] + y.multi.alt[,f.IX]),2)
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#[1,] 0.62 0.57 0.38 0.71 0.28 0.38 0.44 0.41
#[2,] 0.37 0.51 0.48 0.77 0.36 0.35 0.44 0.42
#[3,] 0.53 0.53 0.42 0.64 0.32 0.42 0.66 0.50
#[4,] 0.36 0.55 0.36 0.68 0.35 0.33 0.60 0.57 








##### Use another data set's pi and sd in a different way

## binomial data set <- use multinomial pi and sd

llk.f = loglik_conv(fitted.g.multi.alt[[j]], zdat.bin.alt[3,ind], zdat.bin.alt[4,ind], "+")
matrix_lik = t(compdens_conv(fitted.g.multi.alt[[j]], zdat.bin.alt[3,ind], zdat.bin.alt[4,ind]))
llk.l = sum(log(matrix_lik[,1]))

(llk.f - llk.l)/(2^j)
#[1] -1.202417

## multinomial data set <- use binomial pi and sd

llk.f = loglik_conv(fitted.g.bin.alt[[j]], zdat.multi.alt[3,ind], zdat.multi.alt[4,ind], "+")
matrix_lik = t(compdens_conv(fitted.g.bin.alt[[j]], zdat.multi.alt[3,ind], zdat.multi.alt[4,ind]))
llk.l = sum(log(matrix_lik[,1]))

(llk.f - llk.l)/(2^j)
# 0



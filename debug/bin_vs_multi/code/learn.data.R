



setwd("/mnt/lustre/home/shim/multiscale_analysis")

library("multiseq")
library("ashr")

multiscale.analysis.repodir <- scan(".multiscale_analysis.repodir.txt", what=character())
#source(paste0(multiscale.analysis.repodir, "/src/R/my.utils.R"))

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


locus.start + 2^16
#[1] 25755649
locus.start + 2^16 + 2^15
#[1] 25788417
locus.start + 2^16 + 2^15 + 2^14
#[1] 25804801
peak.start
#[1] 25799098
peak.end
#[1] 25799441

ix17 = locus.start:locus.end
myix16 = (2^16+1):2^17
ix16 = ix17[myix16]
min(ix16)
max(ix16)
peak.start
peak.end


myix15 = 2^16 + (2^15+1):2^16
ix15 = ix17[myix15]
min(ix15)
max(ix15)
peak.start
peak.end


myix14 = 2^16 + 2^15 + 1:(2^14-1)
ix14 = ix17[myix14]
min(ix14)
max(ix14)
peak.start
peak.end


myix13 = 2^16 + 2^15 + 2^13 + 1:(2^13-1)
ix13 = ix17[myix13]
min(ix13)
max(ix13)
peak.start
peak.end


myix12 = 2^16 + 2^15 + 2^13 + 1:(2^12-1)
ix12 = ix17[myix12]
min(ix12)
max(ix12)
peak.start
peak.end




bin.alt.16 = multiseq(x = sim.data.bin.alt[[1]][,myix16], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

bin.null.16 = multiseq(x = sim.data.bin.null[[1]][,myix16], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

multi.alt.16 = multiseq(x = sim.data.multi.alt[[1]][,myix16], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

multi.null.16 = multiseq(x = sim.data.multi.null[[1]][,myix16], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)


bin.alt.16$logLR$value
#[1] 0.1716326
bin.null.16$logLR$value
#[1] 0
multi.alt.16$logLR$value
#[1] 8.32672
multi.null.16$logLR$value
#[1] 0.1088229






bin.alt.15 = multiseq(x = sim.data.bin.alt[[1]][,myix15], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

bin.null.15 = multiseq(x = sim.data.bin.null[[1]][,myix15], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

multi.alt.15 = multiseq(x = sim.data.multi.alt[[1]][,myix15], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

multi.null.15 = multiseq(x = sim.data.multi.null[[1]][,myix15], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)


bin.alt.15$logLR$value
#[1] 0.6159208
bin.null.15$logLR$value
#[1] 0
multi.alt.15$logLR$value
#[1] 10.18559
multi.null.15$logLR$value
#[1] 0.2385997







bin.alt.14 = multiseq(x = sim.data.bin.alt[[1]][,myix14], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

bin.null.14 = multiseq(x = sim.data.bin.null[[1]][,myix14], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

multi.alt.14 = multiseq(x = sim.data.multi.alt[[1]][,myix14], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

multi.null.14 = multiseq(x = sim.data.multi.null[[1]][,myix14], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)


bin.alt.14$logLR$value 
#[1] 2.511342 
bin.null.14$logLR$value
#[1] 0 
multi.alt.14$logLR$value
#[1] 24.36016
multi.null.14$logLR$value
#[1] 0.5599849






bin.alt.13 = multiseq(x = sim.data.bin.alt[[1]][,myix13], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

bin.null.13 = multiseq(x = sim.data.bin.null[[1]][,myix13], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

multi.alt.13 = multiseq(x = sim.data.multi.alt[[1]][,myix13], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

multi.null.13 = multiseq(x = sim.data.multi.null[[1]][,myix13], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)


bin.alt.13$logLR$value 
#[1] 3.927957
bin.null.13$logLR$value
#[1] 0
multi.alt.13$logLR$value
#[1] 28.54315
multi.null.13$logLR$value
#[1] 0.06399062







bin.alt.12 = multiseq(x = sim.data.bin.alt[[1]][,myix12], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

bin.null.12 = multiseq(x = sim.data.bin.null[[1]][,myix12], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

multi.alt.12 = multiseq(x = sim.data.multi.alt[[1]][,myix12], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)

multi.null.12 = multiseq(x = sim.data.multi.null[[1]][,myix12], g = c(0, 0, 1, 1), read.depth = read.depth, prior = "uniform", onlylogLR = TRUE)


bin.alt.12$logLR$value 
#[1] 10.05691
bin.null.12$logLR$value
#[1] 0
multi.alt.12$logLR$value
#[1] 29.66325
multi.null.12$logLR$value
#[1] 0.3078389




#### compare myix15 and myix12

round(bin.null.15$logLR$scales,3)
round(bin.alt.15$logLR$scales,3)
round(multi.null.15$logLR$scales,3)
round(multi.alt.15$logLR$scales,3)

round(bin.null.12$logLR$scales,3)
round(bin.alt.12$logLR$scales,3)
round(multi.null.12$logLR$scales,3)
round(multi.alt.12$logLR$scales,3)

> round(bin.null.15$logLR$scales,3)
 [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
> round(bin.alt.15$logLR$scales,3)
 [1] 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.164 0.000
[13] 0.331 0.121 0.000 0.000
> round(multi.null.15$logLR$scales,3)
 [1] 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.239 0.000 0.000
[13] 0.000 0.000 0.000 0.000
> round(multi.alt.15$logLR$scales,3)
 [1] 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.517 2.515 1.266
[13] 2.978 2.869 0.000 0.041
> round(bin.null.12$logLR$scales,3)
 [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0
> round(bin.alt.12$logLR$scales,3)
 [1] 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.266 0.558 3.642 0.527
[13] 0.000 5.064
> round(multi.null.12$logLR$scales,3)
 [1] 0.000 0.000 0.000 0.000 0.000 0.122 0.000 0.000 0.000 0.000 0.000 0.040
[13] 0.025 0.120
> round(multi.alt.12$logLR$scales,3)
 [1]  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  2.879
[11]  5.940  1.358  1.065 18.421

round(multi.alt.12$fitted.g[[14]]$pi,2)
round(multi.alt.12$fitted.g[[14]]$sd,2)
round(bin.alt.12$fitted.g[[14]]$pi,2)
round(bin.alt.12$fitted.g[[14]]$sd,2)

round(multi.null.12$fitted.g[[14]]$pi,2)
round(multi.null.12$fitted.g[[14]]$sd,2)
round(bin.null.12$fitted.g[[14]]$pi,2)
round(bin.null.12$fitted.g[[14]]$sd,2)


round(multi.alt.12$fitted.g[[11]]$pi,2)
round(multi.alt.12$fitted.g[[11]]$sd,2)
round(bin.alt.12$fitted.g[[11]]$pi,2)
round(bin.alt.12$fitted.g[[11]]$sd,2)

round(multi.null.12$fitted.g[[11]]$pi,2)
round(multi.null.12$fitted.g[[11]]$sd,2)
round(bin.null.12$fitted.g[[11]]$pi,2)
round(bin.null.12$fitted.g[[11]]$sd,2)


> round(multi.alt.12$fitted.g[[14]]$pi,2)
[1] 0 0 0 0 0 0 0 1 0
> round(multi.alt.12$fitted.g[[14]]$sd,2)
[1] 0.00 0.01 0.03 0.05 0.11 0.21 0.42 0.85 1.69
> round(bin.alt.12$fitted.g[[14]]$pi,2)
[1] 0 0 0 0 0 0 0 1 0
> round(bin.alt.12$fitted.g[[14]]$sd,2)
[1] 0.00 0.01 0.02 0.03 0.06 0.13 0.26 0.51 1.02

> round(multi.null.12$fitted.g[[14]]$pi,2)
[1] 0 0 0 0 0 1 0
> round(multi.null.12$fitted.g[[14]]$sd,2)
[1] 0.00 0.01 0.02 0.04 0.09 0.18 0.36
> round(bin.null.12$fitted.g[[14]]$pi,2)
[1] 1 0 0 0 0
> round(bin.null.12$fitted.g[[14]]$sd,2)
[1] 0.00 0.01 0.02 0.05 0.09


> round(multi.alt.12$fitted.g[[11]]$pi,2)
 [1] 0.48 0.00 0.00 0.00 0.00 0.00 0.00 0.52 0.00 0.00
> round(multi.alt.12$fitted.g[[11]]$sd,2)
 [1] 0.00 0.03 0.06 0.12 0.24 0.48 0.95 1.90 3.80 7.61
> round(bin.alt.12$fitted.g[[11]]$pi,2)
[1] 0 0 0 0 0 0 1 0 0
> round(bin.alt.12$fitted.g[[11]]$sd,2)
[1] 0.00 0.03 0.06 0.12 0.23 0.47 0.93 1.86 3.72

> round(multi.null.12$fitted.g[[11]]$pi,2)
[1] 1 0 0 0 0 0 0 0
> round(multi.null.12$fitted.g[[11]]$sd,2)
[1] 0.00 0.02 0.05 0.10 0.19 0.38 0.76 1.53
> round(bin.null.12$fitted.g[[11]]$pi,2)
[1] 1 0 0 0 0 0 0 0
> round(bin.null.12$fitted.g[[11]]$sd,2)
[1] 0.00 0.03 0.05 0.11 0.21 0.43 0.85 1.71




round(multi.alt.15$fitted.g[[14]]$pi,2)
round(multi.alt.15$fitted.g[[14]]$sd,2)
round(bin.alt.15$fitted.g[[14]]$pi,2)
round(bin.alt.15$fitted.g[[14]]$sd,2)

round(multi.null.15$fitted.g[[14]]$pi,2)
round(multi.null.15$fitted.g[[14]]$sd,2)
round(bin.null.15$fitted.g[[14]]$pi,2)
round(bin.null.15$fitted.g[[14]]$sd,2)



round(multi.alt.15$fitted.g[[13]]$pi,2)
round(multi.alt.15$fitted.g[[13]]$sd,2)
round(bin.alt.15$fitted.g[[13]]$pi,2)
round(bin.alt.15$fitted.g[[13]]$sd,2)

round(multi.null.15$fitted.g[[13]]$pi,2)
round(multi.null.15$fitted.g[[13]]$sd,2)
round(bin.null.15$fitted.g[[13]]$pi,2)
round(bin.null.15$fitted.g[[13]]$sd,2)


round(multi.alt.15$fitted.g[[11]]$pi,2)
round(multi.alt.15$fitted.g[[11]]$sd,2)
round(bin.alt.15$fitted.g[[11]]$pi,2)
round(bin.alt.15$fitted.g[[11]]$sd,2)

round(multi.null.15$fitted.g[[11]]$pi,2)
round(multi.null.15$fitted.g[[11]]$sd,2)
round(bin.null.15$fitted.g[[11]]$pi,2)
round(bin.null.15$fitted.g[[11]]$sd,2)



> round(multi.alt.15$fitted.g[[14]]$pi,2)
[1] 0 0 0 0 0 0 1 0 0
> round(multi.alt.15$fitted.g[[14]]$sd,2)
[1] 0.00 0.01 0.03 0.06 0.11 0.22 0.45 0.90 1.80
> round(bin.alt.15$fitted.g[[14]]$pi,2)
[1] 0.00 0.00 0.00 0.00 0.42 0.58 0.00 0.00
> round(bin.alt.15$fitted.g[[14]]$sd,2)
[1] 0.00 0.01 0.03 0.05 0.10 0.20 0.40 0.81
> 
> round(multi.null.15$fitted.g[[14]]$pi,2)
[1] 1 0 0 0 0 0 0 0
> round(multi.null.15$fitted.g[[14]]$sd,2)
[1] 0.00 0.01 0.02 0.05 0.10 0.19 0.39 0.77
> round(bin.null.15$fitted.g[[14]]$pi,2)
[1] 1 0 0 0 0 0 0 0
> round(bin.null.15$fitted.g[[14]]$sd,2)
[1] 0.00 0.01 0.02 0.03 0.07 0.13 0.27 0.54


> round(multi.alt.15$fitted.g[[13]]$pi,2)
 [1] 0.15 0.00 0.00 0.00 0.00 0.00 0.53 0.32 0.00 0.00
> round(multi.alt.15$fitted.g[[13]]$sd,2)
 [1] 0.00 0.01 0.02 0.05 0.10 0.20 0.39 0.78 1.56 3.12
> round(bin.alt.15$fitted.g[[13]]$pi,2)
[1] 0 0 0 0 0 1 0 0 0
> round(bin.alt.15$fitted.g[[13]]$sd,2)
[1] 0.00 0.01 0.03 0.05 0.11 0.21 0.43 0.85 1.70
>

> round(multi.null.15$fitted.g[[13]]$pi,2)
[1] 1 0 0 0 0 0 0 0 0
> round(multi.null.15$fitted.g[[13]]$sd,2)
[1] 0.00 0.01 0.02 0.05 0.09 0.19 0.38 0.76 1.52
> round(bin.null.15$fitted.g[[13]]$pi,2)
[1] 1 0 0 0 0 0 0 0
> round(bin.null.15$fitted.g[[13]]$sd,2)
[1] 0.00 0.01 0.02 0.05 0.09 0.19 0.38 0.76


> round(multi.alt.15$fitted.g[[11]]$pi,2)
 [1] 0.93 0.00 0.00 0.00 0.00 0.00 0.00 0.07 0.00 0.00
> round(multi.alt.15$fitted.g[[11]]$sd,2)
 [1] 0.00 0.03 0.06 0.12 0.24 0.48 0.95 1.90 3.80 7.61
> round(bin.alt.15$fitted.g[[11]]$pi,2)
[1] 0.94 0.00 0.00 0.00 0.00 0.00 0.06 0.00 0.00
> round(bin.alt.15$fitted.g[[11]]$sd,2)
[1] 0.00 0.03 0.06 0.12 0.25 0.49 0.99 1.98 3.95


> round(multi.null.15$fitted.g[[11]]$pi,2)
[1] 1 0 0 0 0 0 0 0 0
> round(multi.null.15$fitted.g[[11]]$sd,2)
[1] 0.00 0.02 0.05 0.09 0.18 0.37 0.74 1.47 2.94
> round(bin.null.15$fitted.g[[11]]$pi,2)
[1] 1 0 0 0 0 0 0 0
> round(bin.null.15$fitted.g[[11]]$sd,2)
[1] 0.00 0.04 0.07 0.15 0.29 0.59 1.18 2.36
 








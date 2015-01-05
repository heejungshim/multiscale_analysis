#!/usr/bin/env Rscript

## Aim : This file contains Rscripts to try different effect size estimates schemes (for simulations).
## To obtain much smoother effect size, we compute mean curve for each group (sig0 and sig1), concatenate them into sig.all (sig.all = c(sig0, sig1); motivation is to allow two mean curves to share the same shrinkage parameters), smooth sig.all using BAYES.THR in wavethresh package with default option, and finally for positions with difference between two smooth mean curves  less then equal to cut.thresh/70, make their difference zero (sig1.smooth[position] = sig0.smooth[position]. Motivation of the last step is that if difference less than equal to cut.thresh/70 roughly corresponds to two groups have difference in sum of read counts less than cut.thresh/2. 
##
## 
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+


chr.list = c( 17 , 17 , 12 , 2, 1,  2,  7,  8, 19,  4,  6,  7,  9, 19)
site.list = c( 570 , 570 , 171, 1617 ,166,  106,  762,  965, 3142, 1258, 1240,  738, 1761, 1668)
genoIX.list =c(11 ,11 , 19, 12, 21, 10,  4,  7,  8, 11, 15,  9, 14,  8)
case.list = 1:14
dir.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/"

ss = 1

#####################
# get a data from Roger's code
#####################
chr = chr.list[ss]
site = site.list[ss]
genoIX = genoIX.list[ss]

path = paste("/mnt/lustre/home/shim/wavelets/data/DNase/region_01_sel_step1/chr", chr, ".loc", sep="")
loc_dat = read.table(path, as.is = TRUE)
st.posi = as.numeric(loc_dat[site.list[ss],2])
en.posi = as.numeric(loc_dat[site.list[ss],3])



## DNase data
path = paste("/mnt/lustre/home/shim/wavelets/data/DNase/sel_data_01_step1/DNase.", chr, ".", site, ".txt", sep="")
DNase_in = read.table(path)
DNase.Roger = DNase_in


## Mappability 
dir_map_path = paste("/mnt/lustre/home/shim/wavelets/data/DNase/sel_data_map_01_step1/DNase.", chr, sep="")
map_path  = paste(dir_map_path, site, "txt", sep=".")
dat_map = read.table(map_path)	
#dim(dat_map)
# 1 2050
map.Roger = dat_map


## genotype
geno.path = paste0("/mnt/lustre/home/shim/wavelets/data/DNase/geno_01_step1/geno_maf/chr", chr, ".", site, ".geno")
genoF = read.table(geno.path, as.is = TRUE)
genoD = genoF[genoIX, 4:73]
genoR = as.numeric(round(genoD))
#table(genoR)
# 0 1 2
# 28 30 12

################################
# combine two strands and
# take into account mappability
################################

numBPs = 1025
numINDs = 70
map = rep(0, numBPs*2)
wh = which(map.Roger[1,] == 1)
map[wh] = rep(1, length(wh))

t_dat = matrix(data = 0, nr = numINDs, nc = numBPs*2)
t_dat[,wh] = as.matrix(DNase.Roger[,wh])

phenoD = t_dat[,1:(numBPs-1)] + t_dat[,(numBPs+1):(numBPs+numBPs-1)]
#com_map = map[1:numBPs] + map[(numBPs+1):(numBPs+numBPs)]
#wh2 = which(com_map > 0)
#av_dat = t(t(com_dat[,wh2])/com_map[wh2])
#com_dat[,wh2] = av_dat


#######################
# alternative signal
#######################


##############################
# Alternative signal
##############################
wh0 = which(genoR == 0)
wh1 = which(genoR == 1)
wh2 = which(genoR == 2)
length(wh0)
# 28
length(wh1)
# 30
length(wh2)
# 12



# Take a mean profile
sig0 = apply(phenoD[wh0,], 2, mean)
sig1 = apply(phenoD[wh1,], 2, mean)


# handle 0 count
#wh.zero = which(sig0 == 0)
#if(length(wh.zero) > 0){ 
#    sig0[wh.zero] = 1/70
#}
#wh.zero = which(sig1 == 0)
#if(length(wh.zero) > 0){ 
#    sig1[wh.zero] = 1/70
#}

## denoise mean curves using wavethresh
library(wavethresh)

## combine two signals
sig.all = c(sig0, sig1)
sig.all.smooth = BAYES.THR(sig.all)

sig0.smooth = sig.all.smooth[1:1024]
sig1.smooth = sig.all.smooth[1025:2048]

for(val in c(2,6,10)){
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


## get smooth signal as a ratio
smooth.ratio = sig1.smooth/sig0.smooth
output = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/data/"
write.table(smooth.ratio, file = paste0(output, "smooth.ratio.", val), col.names = FALSE, row.names = FALSE, quote=FALSE)




##########################################
# Make a figure to see how they look like
##########################################

output = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/fig/effectsize/"
png(paste0(output, "FootprintSignal", val ,".png"), height = 4, width = 7, units="in", res=300)


numBPs = 1024
map = rep(0, numBPs)
wh = which((dat_map[1,1:numBPs] == 1) | (dat_map[1,(1+1025):(numBPs+1025)] == 1))
map[wh] = rep(1, length(wh))
xmin = st.posi
xmax = en.posi -1 
xval = xmin:xmax
xval_mapp = xval[which(map == 1)]


nf <- layout(matrix(1:4,4,1,byrow = TRUE))

#### Group 0
ymax = max(sig0, sig0.smooth, sig1, sig1.smooth)

par(mar=c(3,3,1,1))
plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(0, ymax),axes= FALSE, main ="Major Homozygotes")
axis(1, at=xval_mapp, labels=xval_mapp)
lines(xval, sig0, type = "l", col = "orange")
lines(xval, sig0.smooth, type = "l", col = "red")
axis(2, font=2)
box()


#### Group 1
ymax = max(sig0, sig0.smooth, sig1, sig1.smooth)

par(mar=c(3,3,1,1))
plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(0, ymax),axes= FALSE, main ="Heterozygotes")
axis(1, at=xval_mapp, labels=xval_mapp)
lines(xval, sig1, type = "l", col = "skyblue")
lines(xval, sig1.smooth, type = "l", col = "blue")
axis(2, font=2)
box()


#### smooth signal
ymax = max(sig0, sig0.smooth, sig1, sig1.smooth)

par(mar=c(3,3,1,1))
plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(0, ymax),axes= FALSE, main ="")
axis(1, at=xval_mapp, labels=xval_mapp)
lines(xval, sig0.smooth, type = "l", col = "red")
lines(xval, sig1.smooth, type = "l", col = "blue")
axis(2, font=2)
box()

#### difference
ymax = max(sig0.smooth - sig1.smooth)
ymin = min(sig0.smooth - sig1.smooth)

par(mar=c(3,3,1,1))
plot(1,1,type="n", xlab = "position", ylab = "",xlim = c(xmin, xmax), ylim=c(ymin, ymax),axes= FALSE, main ="Major Homozygotes - Heterozygotes")
axis(1, at=xval_mapp, labels=xval_mapp)
lines(xval, sig0.smooth - sig1.smooth, type = "l", col = "darkgreen")
abline(h = 0)
axis(2, font=2)
box()



dev.off()


}





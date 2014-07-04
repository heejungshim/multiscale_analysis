

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


#library("zoo")

#sliding.win.smooth = function(signal, window.len = 11){
#
#    TS = zoo(signal)
#    out.signal = rollapply(TS, width = window.len, by= 1, FUN = mean, align = "left")
#    resid.len = trunc(window.len/2)
#    res.signal = rep(NA, length(signal))
#    res.signal[1:resid.len] = out.signal[1]
#    res.signal[(resid.len+1):(length(signal)-resid.len)] = out.signal
#    res.signal[(length(signal)-resid.len+1):length(signal)] = out.signal[length(out.signal)]
#    return(res.signal)
#}



##############################
# Alternative signal
##############################
wh0 = which(genoR == 0)
wh1 = which(genoR == 1)
wh2 = which(genoR == 2)
length(wh0)
# 26
length(wh1)
# 27
length(wh2)
# 17



# Take a mean profile
sig0 = apply(phenoD[wh0,], 2, mean)
sig1 = apply(phenoD[wh1,], 2, mean)


# handle 0 count
wh.zero = which(sig0 == 0)
if(length(wh.zero) > 0){ 
    sig0[wh.zero] = 1/70
}
wh.zero = which(sig1 == 0)
if(length(wh.zero) > 0){ 
    sig1[wh.zero] = 1/70
}

# denoise mean curves using wavethresh
library(wavethresh)

w.sig0 = wd(sig0)
w.sig0.smooth = threshold(w.sig0)
sig0.smooth = wr(w.sig0.smooth)

w.sig1 = wd(sig1)
w.sig1.smooth = threshold(w.sig1)
sig1.smooth = wr(w.sig1.smooth)


smooth.ratio = sig1.smooth/sig0.smooth


#path.output = paste0(dir.path, "data", "/raw.dat")
#write.table(phenoD, path.output, row.names = FALSE, col.names = FALSE, quote = FALSE)

#path.output = paste0(dir.path, "data", "/smooth.ratio")
#write.table(smooth.ratio, path.output, row.names = FALSE, col.names = FALSE, quote = FALSE)


# 70 individual
#genoR.f = c(rep(0, 40), rep(1, 30))

#path.output = paste0(dir.path, "data", "/geno70.dat")
#write.table(genoR.f, path.output, row.names = FALSE, col.names = FALSE, quote = FALSE)

#path.output = paste0(dir.path, "data", "/geno70.txt")
#cat("rs A T ", file = path.output)
#cat(genoR.f, file = path.output, append = TRUE)


# 30 individual
#genoR.f = c(rep(0, 15), rep(1, 15))

#path.output = paste0(dir.path, "data", "/geno30.dat")
#write.table(genoR.f, path.output, row.names = FALSE, col.names = FALSE, quote = FALSE)

#path.output = paste0(dir.path, "data", "/geno30.txt")
#cat("rs A T ", file = path.output)
#cat(genoR.f, file = path.output, append = TRUE)



# 30 individual
#genoR.f = c(rep(0, 15), rep(1, 15))

#path.output = paste0(dir.path, "data", "/geno30.dat")
#write.table(genoR.f, path.output, row.names = FALSE, col.names = FALSE, quote = FALSE)

#path.output = paste0(dir.path, "data", "/geno30.txt")
#cat("rs A T ", file = path.output)
#cat(genoR.f, file = path.output, append = TRUE)




# 10 individual
#genoR.f = c(rep(0, 5), rep(1, 5))

#path.output = paste0(dir.path, "data", "/geno10.dat")
#write.table(genoR.f, path.output, row.names = FALSE, col.names = FALSE, quote = FALSE)

#path.output = paste0(dir.path, "data", "/geno10.txt")
#cat("rs A T ", file = path.output)
#cat(genoR.f, file = path.output, append = TRUE)







##########################################
# Make a figure to see how they look like
##########################################

pdf(paste0("footprint_signal",".pdf"), height = 4, width = 7)


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








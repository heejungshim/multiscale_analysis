#!/usr/bin/env Rscript

## Aim : This file contains script to make figures to compare performance of WaveQTL and multiseq using their AuROC over different sample sizes and library read depths. 
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+


###############################
## Read AuROC.wave and AuROC.ms (computed using compute.AUROC.R)
## 5 by 4
## 5 (70, 30, 10, 6, 4)
## 4 (a half, full, 2full, 4full
###############################

input.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/res/AuROC.ms"
AuROC.ms = read.table(file = input.path)[,-1]

input.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/res/AuROC.wave"
AuROC.wave = read.table(file = input.path)[,-1]


############################
## Different sample sizes
############################

pdf("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/code_fig/RECOMB_2014/fig/AuROC_differentSampleSizes.pdf")

plot(0,0, xlim = c(0.5, 1), ylim = c(0.5, 1), xlab = "WaveQTL AuROC", ylab = "multiseq AuROC")
lines(c(-2,2), c(-2,2))
points(AuROC.wave[,2], AuROC.ms[,2], col = c("red", "blue", "purple", "brown", "darkgreen"), pch = 16)
legend(0.5,1, c(70, 30, 10, 6, 4), col = c("red", "blue", "purple", "brown", "darkgreen"), pch = 16, text.col = "black",merge = FALSE, bg = "white")

dev.off()



############################
## Different library read depths 
############################

pdf("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/code_fig/RECOMB_2014/fig/AuROC_differentLibraryReadDepths.pdf")

plot(0,0, xlim = c(0.5, 1), ylim = c(0.5, 1), xlab = "WaveQTL AuROC", ylab = "multiseq AuROC")
lines(c(-2,2), c(-2,2))
points(AuROC.wave[,2], AuROC.ms[,2], col = c("red", "blue", "purple", "gold", "green"), pch = 16)
points(AuROC.wave[1:2,1], AuROC.ms[1:2,1], col = c("red", "blue"), pch = 15)
points(AuROC.wave[3:5,3], AuROC.ms[3:5,3], col = c("purple", "gold", "green"), pch = 17)
points(AuROC.wave[3:5,4], AuROC.ms[3:5,4], col = c("purple", "gold", "green"), pch = 18)
legend(0.5,1, c("70", "30", "10", "6", "4", "50%", "100%", "200%", "400%"), col = c("red", "blue", "purple", "gold", "green", "black", "black", "black", "black"), pch = c(rep(16,5), 15, 16, 17, 18), text.col = "black",merge = FALSE, bg = "white")



dev.off()





############################
## Different sample sizes
############################

pdf("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/code_fig/RECOMB_2014/fig/AuROC_differentSampleSizes_v2.pdf")

d100 = (AuROC.ms[,2] - AuROC.wave[,2])/AuROC.wave[,2]*100 
d50 = (AuROC.ms[1:2,1] - AuROC.wave[1:2,1])/AuROC.wave[1:2,1]*100
d200 = (AuROC.ms[3:5,3] - AuROC.wave[3:5,3])/AuROC.wave[3:5,3]*100
d400 = (AuROC.ms[3:5,4] - AuROC.wave[3:5,4])/AuROC.wave[3:5,4]*100

ymax = max(d100, d50, d200, d400)
ymin = min(d100, d50, d200, d400)


plot(0,0, xlim = c(0.5, 1), ylim = c(-2, 13), xlab = "WaveQTL AuROC", ylab = "fractional change in AuROC")
abline(h=0)
points(AuROC.wave[,2], d100, col = c("red", "blue", "purple", "gold", "green"), pch = 16)
points(AuROC.wave[1:2,1], d50, col = c("red", "blue"), pch = 15)
points(AuROC.wave[3:5,3], d200, col = c("purple", "gold", "green"), pch = 17)
points(AuROC.wave[3:5,4], d400, col = c("purple", "gold", "green"), pch = 18)
legend(0.5,13, c("70", "30", "10", "6", "4", "50%", "100%", "200%", "400%"), col = c("red", "blue", "purple", "gold", "green", "black", "black", "black", "black"), pch = c(rep(16,5), 15, 16, 17, 18), text.col = "black",merge = FALSE, bg = "white")




dev.off()



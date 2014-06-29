


#output.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.Robj"
#output.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.overS.Robj"
#output.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.overB.Robj"
#output.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.overS.RD.Robj"



load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.overS.Robj")

fpr.wave = fpr.wave.list
tpr.wave = tpr.wave.list
fpr.ms = fpr.ms.list
tpr.ms = tpr.ms.list



load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.overS.RD.Robj")

fpr.wave.RD = fpr.wave.list
tpr.wave.RD = tpr.wave.list
fpr.ms.RD = fpr.ms.list
tpr.ms.RD = tpr.ms.list



pdf("logLR_ROCoverS_70.pdf")

xmax = 1
ymax = 1

plot(0,0, xlim = c(0, xmax), ylim = c(0, ymax), xlab = "FPR (1 - Specificity)", ylab ="TPR (Sensitivity)", type="n")
points(c(0, fpr.ms[[1]]), c(0, tpr.ms[[1]]), type="l", col ="red")
points(c(0, fpr.wave[[1]]), c(0, tpr.wave[[1]]), type="l", col ="blue")
#points(c(0, fpr.ms[[2]]), c(0, tpr.ms[[2]]), type="l", col ="red")
#points(c(0, fpr.wave[[2]]), c(0, tpr.wave[[2]]), type="l", col ="blue")
#points(c(0, fpr.ms[[3]]), c(0, tpr.ms[[3]]), type="l", col ="red")
#points(c(0, fpr.wave[[3]]), c(0, tpr.wave[[3]]), type="l", col ="blue")


points(c(0, fpr.ms.RD[[1]]), c(0, tpr.ms.RD[[1]]), type="l", col ="red", lty="dashed")
points(c(0, fpr.wave.RD[[1]]), c(0, tpr.wave.RD[[1]]), type="l", col ="blue", lty="dashed")
#points(c(0, fpr.ms.RD[[2]]), c(0, tpr.ms.RD[[2]]), type="l", col ="red", lty="dashed")
#points(c(0, fpr.wave.RD[[2]]), c(0, tpr.wave.RD[[2]]), type="l", col ="blue", lty="dashed")
#points(c(0, fpr.ms.RD[[3]]), c(0, tpr.ms.RD[[3]]), type="l", col ="red", lty="dashed")
#points(c(0, fpr.wave.RD[[3]]), c(0, tpr.wave.RD[[3]]), type="l", col ="blue", lty="dashed")


legend(0.6, 0.2, c("multiscale (full)", "wave (full)", "multiscale (a half)", "wave (a half)") , col = c("red", "blue", "red", "blue"), lty = c(rep("solid", 2), rep("dashed", 2)), text.col = "black",merge = FALSE, bg = "white")

dev.off()






pdf("logLR_ROCoverS_30.pdf")

xmax = 1
ymax = 1

plot(0,0, xlim = c(0, xmax), ylim = c(0, ymax), xlab = "FPR (1 - Specificity)", ylab ="TPR (Sensitivity)", type="n")
#points(c(0, fpr.ms[[1]]), c(0, tpr.ms[[1]]), type="l", col ="red")
#points(c(0, fpr.wave[[1]]), c(0, tpr.wave[[1]]), type="l", col ="blue")
points(c(0, fpr.ms[[2]]), c(0, tpr.ms[[2]]), type="l", col ="red")
points(c(0, fpr.wave[[2]]), c(0, tpr.wave[[2]]), type="l", col ="blue")
#points(c(0, fpr.ms[[3]]), c(0, tpr.ms[[3]]), type="l", col ="red")
#points(c(0, fpr.wave[[3]]), c(0, tpr.wave[[3]]), type="l", col ="blue")


#points(c(0, fpr.ms.RD[[1]]), c(0, tpr.ms.RD[[1]]), type="l", col ="red", lty="dashed")
#points(c(0, fpr.wave.RD[[1]]), c(0, tpr.wave.RD[[1]]), type="l", col ="blue", lty="dashed")
points(c(0, fpr.ms.RD[[2]]), c(0, tpr.ms.RD[[2]]), type="l", col ="red", lty="dashed")
points(c(0, fpr.wave.RD[[2]]), c(0, tpr.wave.RD[[2]]), type="l", col ="blue", lty="dashed")
#points(c(0, fpr.ms.RD[[3]]), c(0, tpr.ms.RD[[3]]), type="l", col ="red", lty="dashed")
#points(c(0, fpr.wave.RD[[3]]), c(0, tpr.wave.RD[[3]]), type="l", col ="blue", lty="dashed")


legend(0.6, 0.2, c("multiscale (full)", "wave (full)", "multiscale (a half)", "wave (a half)") , col = c("red", "blue", "red", "blue"), lty = c(rep("solid", 2), rep("dashed", 2)), text.col = "black",merge = FALSE, bg = "white")

dev.off()






pdf("logLR_ROCoverS_10.pdf")

xmax = 1
ymax = 1

plot(0,0, xlim = c(0, xmax), ylim = c(0, ymax), xlab = "FPR (1 - Specificity)", ylab ="TPR (Sensitivity)", type="n")
#points(c(0, fpr.ms[[1]]), c(0, tpr.ms[[1]]), type="l", col ="red")
#points(c(0, fpr.wave[[1]]), c(0, tpr.wave[[1]]), type="l", col ="blue")
#points(c(0, fpr.ms[[2]]), c(0, tpr.ms[[2]]), type="l", col ="red")
#points(c(0, fpr.wave[[2]]), c(0, tpr.wave[[2]]), type="l", col ="blue")
points(c(0, fpr.ms[[3]]), c(0, tpr.ms[[3]]), type="l", col ="red")
points(c(0, fpr.wave[[3]]), c(0, tpr.wave[[3]]), type="l", col ="blue")


#points(c(0, fpr.ms.RD[[1]]), c(0, tpr.ms.RD[[1]]), type="l", col ="red", lty="dashed")
#points(c(0, fpr.wave.RD[[1]]), c(0, tpr.wave.RD[[1]]), type="l", col ="blue", lty="dashed")
#points(c(0, fpr.ms.RD[[2]]), c(0, tpr.ms.RD[[2]]), type="l", col ="red", lty="dashed")
#points(c(0, fpr.wave.RD[[2]]), c(0, tpr.wave.RD[[2]]), type="l", col ="blue", lty="dashed")
points(c(0, fpr.ms.RD[[3]]), c(0, tpr.ms.RD[[3]]), type="l", col ="red", lty="dashed")
points(c(0, fpr.wave.RD[[3]]), c(0, tpr.wave.RD[[3]]), type="l", col ="blue", lty="dashed")


legend(0.6, 0.2, c("multiscale (full)", "wave (full)", "multiscale (2 full)", "wave (2 full)") , col = c("red", "blue", "red", "blue"), lty = c(rep("solid", 2), rep("dashed", 2)), text.col = "black",merge = FALSE, bg = "white")

dev.off()








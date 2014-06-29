


#output.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.Robj"
#output.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.overS.Robj"
#output.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.overB.Robj"
#output.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.overS.RD.Robj"



load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.Robj")

fpr.wave = fpr.wave.list
tpr.wave = tpr.wave.list
fpr.ms = fpr.ms.list
tpr.ms = tpr.ms.list


load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.overS.Robj")

fpr.wave.overS = fpr.wave.list
tpr.wave.overS = tpr.wave.list
fpr.ms.overS = fpr.ms.list
tpr.ms.overS = tpr.ms.list



load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.overB.Robj")

fpr.wave.overB = fpr.wave.list
tpr.wave.overB = tpr.wave.list
fpr.ms.overB = fpr.ms.list
tpr.ms.overB = tpr.ms.list


pdf("logLR_ROC_withover.pdf")

xmax = 1
ymax = 1

plot(0,0, xlim = c(0, xmax), ylim = c(0, ymax), xlab = "FPR (1 - Specificity)", ylab ="TPR (Sensitivity)", type="n")
points(c(0, fpr.ms[[1]]), c(0, tpr.ms[[1]]), type="l", col ="blue")
points(c(0, fpr.wave[[1]]), c(0, tpr.wave[[1]]), type="l", col ="skyblue")
points(c(0, fpr.ms[[2]]), c(0, tpr.ms[[2]]), type="l", col ="darkgreen")
points(c(0, fpr.wave[[2]]), c(0, tpr.wave[[2]]), type="l", col ="green")
points(c(0, fpr.ms[[3]]), c(0, tpr.ms[[3]]), type="l", col ="red")
points(c(0, fpr.wave[[3]]), c(0, tpr.wave[[3]]), type="l", col ="orange")

points(c(0, fpr.ms.overS[[1]]), c(0, tpr.ms.overS[[1]]), type="l", col ="blue", lty="twodash")
points(c(0, fpr.wave.overS[[1]]), c(0, tpr.wave.overS[[1]]), type="l", col ="skyblue", lty="twodash")
points(c(0, fpr.ms.overS[[2]]), c(0, tpr.ms.overS[[2]]), type="l", col ="darkgreen", lty="twodash")
points(c(0, fpr.wave.overS[[2]]), c(0, tpr.wave.overS[[2]]), type="l", col ="green", lty="twodash")
points(c(0, fpr.ms.overS[[3]]), c(0, tpr.ms.overS[[3]]), type="l", col ="red", lty="twodash")
points(c(0, fpr.wave.overS[[3]]), c(0, tpr.wave.overS[[3]]), type="l", col ="orange", lty="twodash")


points(c(0, fpr.ms.overB[[1]]), c(0, tpr.ms.overB[[1]]), type="l", col ="blue", lty="dashed")
points(c(0, fpr.wave.overB[[1]]), c(0, tpr.wave.overB[[1]]), type="l", col ="skyblue", lty="dashed")
points(c(0, fpr.ms.overB[[2]]), c(0, tpr.ms.overB[[2]]), type="l", col ="darkgreen", lty="dashed")
points(c(0, fpr.wave.overB[[2]]), c(0, tpr.wave.overB[[2]]), type="l", col ="green", lty="dashed")
points(c(0, fpr.ms.overB[[3]]), c(0, tpr.ms.overB[[3]]), type="l", col ="red", lty="dashed")
points(c(0, fpr.wave.overB[[3]]), c(0, tpr.wave.overB[[3]]), type="l", col ="orange", lty="dashed")


legend(0.6, 0.4, c("multiscale 70", "wave 70", "multiscale 30", "wave 30","multiscale 10", "wave 10", "binomial", "beta-binomial (small)", "beta-binomial (big)") , col = c("blue", "skyblue", "darkgreen", "green", "red", "orange", "black", "black", "black"), lty = c(rep("solid", 7), "twodash", "dashed"), text.col = "black",merge = FALSE, bg = "white")

dev.off()








pdf("logLR_ROC_withoverS.pdf")

xmax = 1
ymax = 1

plot(0,0, xlim = c(0, xmax), ylim = c(0, ymax), xlab = "FPR (1 - Specificity)", ylab ="TPR (Sensitivity)", type="n")
points(c(0, fpr.ms[[1]]), c(0, tpr.ms[[1]]), type="l", col ="blue")
points(c(0, fpr.wave[[1]]), c(0, tpr.wave[[1]]), type="l", col ="skyblue")
points(c(0, fpr.ms[[2]]), c(0, tpr.ms[[2]]), type="l", col ="darkgreen")
points(c(0, fpr.wave[[2]]), c(0, tpr.wave[[2]]), type="l", col ="green")
points(c(0, fpr.ms[[3]]), c(0, tpr.ms[[3]]), type="l", col ="red")
points(c(0, fpr.wave[[3]]), c(0, tpr.wave[[3]]), type="l", col ="orange")

points(c(0, fpr.ms.overS[[1]]), c(0, tpr.ms.overS[[1]]), type="l", col ="blue", lty="twodash")
points(c(0, fpr.wave.overS[[1]]), c(0, tpr.wave.overS[[1]]), type="l", col ="skyblue", lty="twodash")
points(c(0, fpr.ms.overS[[2]]), c(0, tpr.ms.overS[[2]]), type="l", col ="darkgreen", lty="twodash")
points(c(0, fpr.wave.overS[[2]]), c(0, tpr.wave.overS[[2]]), type="l", col ="green", lty="twodash")
points(c(0, fpr.ms.overS[[3]]), c(0, tpr.ms.overS[[3]]), type="l", col ="red", lty="twodash")
points(c(0, fpr.wave.overS[[3]]), c(0, tpr.wave.overS[[3]]), type="l", col ="orange", lty="twodash")


#points(c(0, fpr.ms.overB[[1]]), c(0, tpr.ms.overB[[1]]), type="l", col ="blue", lty="dashed")
#points(c(0, fpr.wave.overB[[1]]), c(0, tpr.wave.overB[[1]]), type="l", col ="skyblue", lty="dashed")
#points(c(0, fpr.ms.overB[[2]]), c(0, tpr.ms.overB[[2]]), type="l", col ="darkgreen", lty="dashed")
#points(c(0, fpr.wave.overB[[2]]), c(0, tpr.wave.overB[[2]]), type="l", col ="green", lty="dashed")
#points(c(0, fpr.ms.overB[[3]]), c(0, tpr.ms.overB[[3]]), type="l", col ="red", lty="dashed")
#points(c(0, fpr.wave.overB[[3]]), c(0, tpr.wave.overB[[3]]), type="l", col ="orange", lty="dashed")


legend(0.6, 0.4, c("multiscale 70", "wave 70", "multiscale 30", "wave 30","multiscale 10", "wave 10", "binomial", "beta-binomial (small)") , col = c("blue", "skyblue", "darkgreen", "green", "red", "orange", "black", "black"), lty = c(rep("solid", 7), "twodash"), text.col = "black",merge = FALSE, bg = "white")

dev.off()








pdf("logLR_ROC_withoverB.pdf")

xmax = 1
ymax = 1

plot(0,0, xlim = c(0, xmax), ylim = c(0, ymax), xlab = "FPR (1 - Specificity)", ylab ="TPR (Sensitivity)", type="n")
points(c(0, fpr.ms[[1]]), c(0, tpr.ms[[1]]), type="l", col ="blue")
points(c(0, fpr.wave[[1]]), c(0, tpr.wave[[1]]), type="l", col ="skyblue")
points(c(0, fpr.ms[[2]]), c(0, tpr.ms[[2]]), type="l", col ="darkgreen")
points(c(0, fpr.wave[[2]]), c(0, tpr.wave[[2]]), type="l", col ="green")
points(c(0, fpr.ms[[3]]), c(0, tpr.ms[[3]]), type="l", col ="red")
points(c(0, fpr.wave[[3]]), c(0, tpr.wave[[3]]), type="l", col ="orange")

#points(c(0, fpr.ms.overS[[1]]), c(0, tpr.ms.overS[[1]]), type="l", col ="blue", lty="twodash")
#points(c(0, fpr.wave.overS[[1]]), c(0, tpr.wave.overS[[1]]), type="l", col ="skyblue", lty="twodash")
#points(c(0, fpr.ms.overS[[2]]), c(0, tpr.ms.overS[[2]]), type="l", col ="darkgreen", lty="twodash")
#points(c(0, fpr.wave.overS[[2]]), c(0, tpr.wave.overS[[2]]), type="l", col ="green", lty="twodash")
#points(c(0, fpr.ms.overS[[3]]), c(0, tpr.ms.overS[[3]]), type="l", col ="red", lty="twodash")
#points(c(0, fpr.wave.overS[[3]]), c(0, tpr.wave.overS[[3]]), type="l", col ="orange", lty="twodash")


points(c(0, fpr.ms.overB[[1]]), c(0, tpr.ms.overB[[1]]), type="l", col ="blue", lty="dashed")
points(c(0, fpr.wave.overB[[1]]), c(0, tpr.wave.overB[[1]]), type="l", col ="skyblue", lty="dashed")
points(c(0, fpr.ms.overB[[2]]), c(0, tpr.ms.overB[[2]]), type="l", col ="darkgreen", lty="dashed")
points(c(0, fpr.wave.overB[[2]]), c(0, tpr.wave.overB[[2]]), type="l", col ="green", lty="dashed")
points(c(0, fpr.ms.overB[[3]]), c(0, tpr.ms.overB[[3]]), type="l", col ="red", lty="dashed")
points(c(0, fpr.wave.overB[[3]]), c(0, tpr.wave.overB[[3]]), type="l", col ="orange", lty="dashed")


legend(0.6, 0.4, c("multiscale 70", "wave 70", "multiscale 30", "wave 30","multiscale 10", "wave 10", "binomial", "beta-binomial (big)") , col = c("blue", "skyblue", "darkgreen", "green", "red", "orange", "black", "black"), lty = c(rep("solid", 7), "dashed"), text.col = "black",merge = FALSE, bg = "white")

dev.off()





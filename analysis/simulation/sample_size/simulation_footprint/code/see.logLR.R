


#case.name = c("fullread.70ind", "fullread.30ind", "fullread.10ind")
#case.name = c("fullread.70ind.over", "fullread.30ind.over", "fullread.10ind.over")

#case.name = c("fullread.70ind.over.2", "fullread.30ind.over.2", "fullread.10ind.over.2")
case.name = c("halfread.70ind.over", "halfread.30ind.over", "2fullread.10ind.over")


#output.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.Robj"
#output.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.overS.Robj"
#output.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.overB.Robj"
output.path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/sum/logLR.f.overS.RD.Robj"




#ROC.file.name = "logLR_ROC.pdf"
#ROC.file.name = "logLR_ROC_overS.pdf"
#ROC.file.name = "logLR_ROC_overB.pdf"
ROC.file.name = "logLR_ROC_overS_RD.pdf"

#hist.file.name = "logLR_hist.pdf"
#hist.file.name = "logLR_hist_overS.pdf"
#hist.file.name = "logLR_hist_overB.pdf"
hist.file.name = "logLR_hist_overS_RD.pdf"



##### ms ######
ms.null = vector("list", length(case.name))
ms.alt = vector("list", length(case.name))


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[1], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[1], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[1]] = logLR.null
ms.alt[[1]] = logLR.alt



load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[2], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[2], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res





sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[2]] = logLR.null
ms.alt[[2]] = logLR.alt




load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/logLR.", case.name[3], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/logLR.", case.name[3], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res




sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

ms.null[[3]] = logLR.null
ms.alt[[3]] = logLR.alt







##### wave ######


wave.null = vector("list", length(case.name))
wave.alt = vector("list", length(case.name))


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[1], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[1], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[1]] = logLR.null
wave.alt[[1]] = logLR.alt


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[2], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[2], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res




sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[2]] = logLR.null
wave.alt[[2]] = logLR.alt



load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/logLR.", case.name[3], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/logLR.", case.name[3], ".Robj"))
logLR.null  = as.numeric(logLR_list)
done.null = done_res



sum(done.alt)
sum(done.null)
min(logLR.alt)
max(logLR.alt)
min(logLR.null)
max(logLR.null)

wave.null[[3]] = logLR.null
wave.alt[[3]] = logLR.alt




#############################
# Make histogram
#############################


pdf(hist.file.name)


par(mfrow = c(4,1))
maxval = max(wave.null[[1]], ms.null[[1]], wave.alt[[1]], ms.alt[[1]])
minval = min(wave.null[[1]], ms.null[[1]], wave.alt[[1]], ms.alt[[1]])

hist(wave.null[[1]], main="null Wavelet (70)", breaks = 47, xlim =c(minval, maxval))
hist(ms.null[[1]], main="null multiscale (70)", breaks = 47, xlim =c(minval, maxval))
hist(wave.alt[[1]], main="alt Wavelet (70)", breaks = 47, xlim =c(minval, maxval))
hist(ms.alt[[1]], main="alt multiscale (70)", breaks = 47, xlim =c(minval, maxval))


par(mfrow = c(4,1))
maxval = max(wave.null[[2]], ms.null[[2]], wave.alt[[2]], ms.alt[[2]])
minval = min(wave.null[[2]], ms.null[[2]], wave.alt[[2]], ms.alt[[2]])

hist(wave.null[[2]], main="null Wavelet (30)", breaks = 47, xlim =c(minval, maxval))
hist(ms.null[[2]], main="null multiscale (30)", breaks = 47, xlim =c(minval, maxval))
hist(wave.alt[[2]], main="alt Wavelet (30)", breaks = 47, xlim =c(minval, maxval))
hist(ms.alt[[2]], main="alt multiscale (30)", breaks = 47, xlim =c(minval, maxval))



par(mfrow = c(4,1))
maxval = max(wave.null[[3]], ms.null[[3]], wave.alt[[3]], ms.alt[[3]])
minval = min(wave.null[[3]], ms.null[[3]], wave.alt[[3]], ms.alt[[3]])

hist(wave.null[[3]], main="null Wavelet (10)", breaks = 47, xlim =c(minval, maxval))
hist(ms.null[[3]], main="null multiscale (10)", breaks = 47, xlim =c(minval, maxval))
hist(wave.alt[[3]], main="alt Wavelet (10)", breaks = 47, xlim =c(minval, maxval))
hist(ms.alt[[3]], main="alt multiscale (10)", breaks = 47, xlim =c(minval, maxval))



dev.off()





 

######## roc curves

fpr.wave.list = vector("list", length(case.name))
tpr.wave.list = vector("list", length(case.name))

fpr.ms.list = vector("list", length(case.name))
tpr.ms.list = vector("list", length(case.name))


## wave
for(cc in 1:length(case.name)){

    #cc = 1
    logLR = as.numeric(c(wave.null[[cc]], wave.alt[[cc]]))
    disc = c(rep(0,500), rep(1,500))
    rnk = order(logLR, decreasing = TRUE)
    p.wave = logLR[rnk]
    d.wave = disc[rnk]

    fdp.wave = NULL
    sig.wave = NULL
    tpr.wave = NULL
    fpr.wave = NULL
    uni.p.wave = unique(p.wave)
    for(i in 1:length(uni.p.wave)){
        wh = which(p.wave >= uni.p.wave[i])
        sig.wave[i] = length(wh)
        fdp.wave[i] = 1 - (sum(d.wave[wh])/length(wh))
        tpr.wave[i] = sum(d.wave[wh])/500
        fpr.wave[i] = (length(wh) - sum(d.wave[wh]))/500
    }

    fpr.wave.list[[cc]] = fpr.wave
    tpr.wave.list[[cc]] = tpr.wave
    
}




## ms
for(cc in 1:length(case.name)){
    
    logLR = as.numeric(c(ms.null[[cc]], ms.alt[[cc]]))
    disc = c(rep(0,500), rep(1,500))
    rnk = order(logLR, decreasing = TRUE)
    p.ms = logLR[rnk]
    d.ms = disc[rnk]

    fdp.ms = NULL
    sig.ms = NULL
    tpr.ms = NULL
    fpr.ms = NULL
    uni.p.ms = unique(p.ms)
    for(i in 1:length(uni.p.ms)){
        wh = which(p.ms >= uni.p.ms[i])
        sig.ms[i] = length(wh)
        fdp.ms[i] = 1 - (sum(d.ms[wh])/length(wh))
        tpr.ms[i] = sum(d.ms[wh])/500
        fpr.ms[i] = (length(wh) - sum(d.ms[wh]))/500
    }

    fpr.ms.list[[cc]] = fpr.ms
    tpr.ms.list[[cc]] = tpr.ms
    
}



save("fpr.wave.list", "tpr.wave.list", "fpr.ms.list", "tpr.ms.list", file = output.path)






pdf(ROC.file.name)



xmax = 1
ymax = 1

plot(0,0, xlim = c(0, xmax), ylim = c(0, ymax), xlab = "FPR (1 - Specificity)", ylab ="TPR (Sensitivity)", type="n")
points(c(0, fpr.ms.list[[1]]), c(0, tpr.ms.list[[1]]), type="l", col ="blue")
points(c(0, fpr.wave.list[[1]]), c(0, tpr.wave.list[[1]]), type="l", col ="skyblue")
points(c(0, fpr.ms.list[[2]]), c(0, tpr.ms.list[[2]]), type="l", col ="darkgreen")
points(c(0, fpr.wave.list[[2]]), c(0, tpr.wave.list[[2]]), type="l", col ="green")
points(c(0, fpr.ms.list[[3]]), c(0, tpr.ms.list[[3]]), type="l", col ="red")
points(c(0, fpr.wave.list[[3]]), c(0, tpr.wave.list[[3]]), type="l", col ="orange")

legend(0.7, 0.3, c("multiscale 70", "wave 70", "multiscale 30", "wave 30","multiscale 10", "wave 10") , col = c("blue", "skyblue", "darkgreen", "green", "red", "orange"), lty = c(1,1), text.col = "black",merge = FALSE, bg = "white")

dev.off()



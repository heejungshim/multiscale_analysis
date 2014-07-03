
numIND = 4

case.name = c(paste0("fullread.", numIND, "ind"), paste0("fullread.", numIND, "ind.over"), paste0("fullread.", numIND, "ind.over.2"))


ROC.file.name = paste0("logLR_ROC_over", numIND, ".pdf")
hist.file.name = paste0("logLR_hist_over", numIND, ".pdf")





##### ms ######
ms.null = vector("list", length(case.name))
ms.alt = vector("list", length(case.name))


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/multiscale/sum/logLR.", case.name[1], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/null/multiscale/sum/logLR.", case.name[1], ".Robj"))
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



load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/multiscale/sum/logLR.", case.name[2], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/null/multiscale/sum/logLR.", case.name[2], ".Robj"))
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




load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/multiscale/sum/logLR.", case.name[3], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/null/multiscale/sum/logLR.", case.name[3], ".Robj"))
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


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/wave/sum/logLR.", case.name[1], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/null/wave/sum/logLR.", case.name[1], ".Robj"))
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


load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/wave/sum/logLR.", case.name[2], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/null/wave/sum/logLR.", case.name[2], ".Robj"))
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



load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/wave/sum/logLR.", case.name[3], ".Robj"))
logLR.alt  = as.numeric(logLR_list)
done.alt = done_res

load(paste0("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/null/wave/sum/logLR.", case.name[3], ".Robj"))
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
hist(wave.null[[1]], main=paste0("null Wavelet : ", numIND), breaks = 47)
hist(ms.null[[1]], main=paste0("null multiscale : ", numIND), breaks = 47)
hist(wave.alt[[1]], main=paste0("alt Wavelet : ", numIND), breaks = 47)
hist(ms.alt[[1]], main=paste0("alt multiscale : ", numIND), breaks = 47)


par(mfrow = c(4,1))
hist(wave.null[[2]], main=paste0("null Wavelet : ", numIND, " OD small"), breaks = 47)
hist(ms.null[[2]], main=paste0("null multiscale : ", numIND, " OD small"), breaks = 47)
hist(wave.alt[[2]], main=paste0("alt Wavelet : ", numIND, " OD small"), breaks = 47)
hist(ms.alt[[2]], main=paste0("alt multiscale : ", numIND, " OD small"), breaks = 47)


par(mfrow = c(4,1))
hist(wave.null[[3]], main=paste0("null Wavelet : ", numIND, " OD big"), breaks = 47)
hist(ms.null[[3]], main=paste0("null multiscale : ", numIND, " OD big"), breaks = 47)
hist(wave.alt[[3]], main=paste0("alt Wavelet : ", numIND, " OD big"), breaks = 47)
hist(ms.alt[[3]], main=paste0("alt multiscale : ", numIND, " OD big"), breaks = 47)



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
    disc = c(rep(0,578), rep(1,578))
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
        tpr.wave[i] = sum(d.wave[wh])/578
        fpr.wave[i] = (length(wh) - sum(d.wave[wh]))/578
    }

    fpr.wave.list[[cc]] = fpr.wave
    tpr.wave.list[[cc]] = tpr.wave
    
}




## ms
for(cc in 1:length(case.name)){
    
    logLR = as.numeric(c(ms.null[[cc]], ms.alt[[cc]]))
    disc = c(rep(0,578), rep(1,578))
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
        tpr.ms[i] = sum(d.ms[wh])/578
        fpr.ms[i] = (length(wh) - sum(d.ms[wh]))/578
    }

    fpr.ms.list[[cc]] = fpr.ms
    tpr.ms.list[[cc]] = tpr.ms
    
}






pdf(ROC.file.name)

xmax = 1
ymax = 1


plot(0,0, xlim = c(0, xmax), ylim = c(0, ymax), xlab = "FPR (1 - Specificity)", ylab ="TPR (Sensitivity)", type="n")
points(c(0, fpr.ms.list[[1]]), c(0, tpr.ms.list[[1]]), type="l", col ="red")
points(c(0, fpr.wave.list[[1]]), c(0, tpr.wave.list[[1]]), type="l", col ="blue")
points(c(0, fpr.ms.list[[2]]), c(0, tpr.ms.list[[2]]), type="l", col ="red", lty="dashed")
points(c(0, fpr.wave.list[[2]]), c(0, tpr.wave.list[[2]]), type="l", col ="blue", lty="dashed")
points(c(0, fpr.ms.list[[3]]), c(0, tpr.ms.list[[3]]), type="l", col ="red", lty="dotdash")
points(c(0, fpr.wave.list[[3]]), c(0, tpr.wave.list[[3]]), type="l", col ="blue", lty="dotdash")

legend(0.6, 0.2, c(paste0("multiscale ", numIND), paste0("Wavelets ", numIND), "beta-binomial small", "beta-binomial big") , col = c("red", "blue", "black", "black"), lty = c("solid", "solid", "dashed", "dotdash" ), text.col = "black",merge = FALSE, bg = "white")

dev.off()



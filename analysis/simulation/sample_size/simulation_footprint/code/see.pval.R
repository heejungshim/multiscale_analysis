


case.name = c("fullread.70ind", "fullread.30ind", "fullread.10ind")

##### ms ######
ms.null = vector("list", length(case.name))
ms.alt = vector("list", length(case.name))


load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/pval.fullread.70ind.Robj")
pval.alt  = as.numeric(pval_list)
done.alt = done_res

load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/pval.fullread.70ind.Robj")
pval.null  = as.numeric(pval_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(pval.alt)
max(pval.alt)
min(pval.null)
max(pval.null)

ms.null[[1]] = pval.null
ms.alt[[1]] = pval.alt


load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/pval.fullread.30ind.Robj")
pval.alt  = as.numeric(pval_list)
done.alt = done_res

load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/pval.fullread.30ind.Robj")
pval.null  = as.numeric(pval_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(pval.alt)
max(pval.alt)
min(pval.null)
max(pval.null)

ms.null[[2]] = pval.null
ms.alt[[2]] = pval.alt


load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/multiscale/sum/pval.fullread.10ind.Robj")
pval.alt  = as.numeric(pval_list)
done.alt = done_res

load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/sum/pval.fullread.10ind.Robj")
pval.null  = as.numeric(pval_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(pval.alt)
max(pval.alt)
min(pval.null)
max(pval.null)

ms.null[[3]] = pval.null
ms.alt[[3]] = pval.alt







##### wave ######


wave.null = vector("list", length(case.name))
wave.alt = vector("list", length(case.name))


load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/pval.fullread.70ind.Robj")
pval.alt  = as.numeric(pval_list)
done.alt = done_res

load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/pval.fullread.70ind.Robj")
pval.null  = as.numeric(pval_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(pval.alt)
max(pval.alt)
min(pval.null)
max(pval.null)

wave.null[[1]] = pval.null
wave.alt[[1]] = pval.alt


load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/pval.fullread.30ind.Robj")
pval.alt  = as.numeric(pval_list)
done.alt = done_res

load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/pval.fullread.30ind.Robj")
pval.null  = as.numeric(pval_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(pval.alt)
max(pval.alt)
min(pval.null)
max(pval.null)

wave.null[[2]] = pval.null
wave.alt[[2]] = pval.alt


load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/sum/pval.fullread.10ind.Robj")
pval.alt  = as.numeric(pval_list)
done.alt = done_res

load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/sum/pval.fullread.10ind.Robj")
pval.null  = as.numeric(pval_list)
done.null = done_res


sum(done.alt)
sum(done.null)
min(pval.alt)
max(pval.alt)
min(pval.null)
max(pval.null)

wave.null[[3]] = pval.null
wave.alt[[3]] = pval.alt







#############################
# Make histogram
#############################

pdf("hist_null_alt_footprint.pdf")

par(mfrow = c(4,1))
hist(wave.null[[1]], main="null Wavelet (70)", breaks = 47)
hist(ms.null[[1]], main="null multiscale (70)", breaks = 47)
hist(wave.alt[[1]], main="alt Wavelet (70)", breaks = 47)
hist(ms.alt[[1]], main="alt multiscale (70)", breaks = 47)


par(mfrow = c(4,1))
hist(wave.null[[2]], main="null Wavelet (30)", breaks = 47)
hist(ms.null[[2]], main="null multiscale (30)", breaks = 47)
hist(wave.alt[[2]], main="alt Wavelet (30)", breaks = 47)
hist(ms.alt[[2]], main="alt multiscale (30)", breaks = 47)



par(mfrow = c(4,1))
hist(wave.null[[3]], main="null Wavelet (10)", breaks = 47)
hist(ms.null[[3]], main="null multiscale (10)", breaks = 47)
hist(wave.alt[[3]], main="alt Wavelet (10)", breaks = 47)
hist(ms.alt[[3]], main="alt multiscale (10)", breaks = 47)



dev.off()






######## roc curves

fpr.wave.list = vector("list", length(case.name))
tpr.wave.list = vector("list", length(case.name))

fpr.ms.list = vector("list", length(case.name))
tpr.ms.list = vector("list", length(case.name))


## wave
for(cc in 1:length(case.name)){
    
    pval = as.numeric(c(wave.null[[cc]], wave.alt[[cc]]))
    disc = c(rep(0,200), rep(1,200))
    rnk = order(pval)
    p.wave = pval[rnk]
    d.wave = disc[rnk]

    fdp.wave = NULL
    sig.wave = NULL
    tpr.wave = NULL
    fpr.wave = NULL
    uni.p.wave = unique(p.wave)
    for(i in 1:length(uni.p.wave)){
        wh = which(p.wave <= uni.p.wave[i])
        sig.wave[i] = length(wh)
        fdp.wave[i] = 1 - (sum(d.wave[wh])/length(wh))
        tpr.wave[i] = sum(d.wave[wh])/200
        fpr.wave[i] = (length(wh) - sum(d.wave[wh]))/200
    }

    fpr.wave.list[[cc]] = fpr.wave
    tpr.wave.list[[cc]] = tpr.wave
    
}




## ms
for(cc in 1:length(case.name)){
    
    pval = as.numeric(c(ms.null[[cc]], ms.alt[[cc]]))
    disc = c(rep(0,200), rep(1,200))
    rnk = order(pval)
    p.ms = pval[rnk]
    d.ms = disc[rnk]

    fdp.ms = NULL
    sig.ms = NULL
    tpr.ms = NULL
    fpr.ms = NULL
    uni.p.ms = unique(p.ms)
    for(i in 1:length(uni.p.ms)){
        wh = which(p.ms <= uni.p.ms[i])
        sig.ms[i] = length(wh)
        fdp.ms[i] = 1 - (sum(d.ms[wh])/length(wh))
        tpr.ms[i] = sum(d.ms[wh])/200
        fpr.ms[i] = (length(wh) - sum(d.ms[wh]))/200
    }

    fpr.ms.list[[cc]] = fpr.ms
    tpr.ms.list[[cc]] = tpr.ms
    
}




pdf("simu_roc_footprint_50.pdf")

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









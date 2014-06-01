


##### Null wave JAR ######
load("/mnt/lustre/home/shim/wavelets/revision/analysis/simulation_578/null/wave/sum/pval.Robj")
pval.null.wave = as.numeric(pval_list)
done.null.wave = done_res

load("/mnt/lustre/home/shim/wavelets/revision/analysis/simulation_578/null/JAR/sum/pval.Robj")
pval.null.JAR = pval_list
done.null.JAR = done_res

sum(done.null.wave)
# 100
sum(done.null.JAR)
# 100



##### alt1 wave JAR ######
load("/mnt/lustre/home/shim/wavelets/revision/analysis/simulation_578/alt/wave/sum_21_80/pval.Robj")
pval.alt.wave = as.numeric(pval_list)
done.alt.wave = done_res


load("/mnt/lustre/home/shim/wavelets/revision/analysis/simulation_578/alt/JAR/sum_21_80/pval.Robj")
pval.alt.JAR = pval_list
done.alt.JAR = done_res

###### multi-scale null and alt


load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/null/multiscale/sum/pval.Robj")
pval.null.ms = as.numeric(pval_list)
done.null.ms = done_res

load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/multiscale/sum/pval.Robj")
pval.alt.ms = as.numeric(pval_list)
done.alt.ms = done_res



pdf("hist_null_alt_578.pdf")

par(mfrow = c(3,1))
hist(pval.null.wave, main="null Wavelet", breaks = 47)
hist(pval.null.JAR, main="null 100-window", breaks = 47)
hist(pval.null.ms, main="null multi-scale", breaks = 47)

par(mfrow = c(3,1))
hist(pval.alt.wave, main="null Wavelet", breaks = 47)
hist(pval.alt.JAR, main="null 100-window", breaks = 47)
hist(pval.alt.ms, main="null multi-scale", breaks = 47)

delIX = c(305, 368, 451, 464)
lp.wave = -log(pval.alt.wave,10)[-delIX]
lp.ms = -log(pval.alt.ms,10)[-delIX]
xymax = max(lp.wave, lp.ms)
plot(pval.alt.wave, pval.alt.ms, xlim=c(0,xymax), ylim=c(0,xymax))

dev.off()


######## FDP

pval.null.ms = pval.null.ms[-delIX]
pval.null.wave = pval.null.wave[-delIX]
pval.null.JAR = pval.null.JAR[-delIX]

pval.alt.ms = pval.alt.ms[-delIX]
pval.alt.wave = pval.alt.wave[-delIX]
pval.alt.JAR = pval.alt.JAR[-delIX]


## for wave
pval.wave = as.numeric(c(pval.null.wave, pval.alt.wave))
disc.wave = c(rep(0,574), rep(1,574))
rnk.wave = order(pval.wave)
p.wave = pval.wave[rnk.wave]
d.wave = disc.wave[rnk.wave]

fdp.wave = NULL
sig.wave = NULL
tpr.wave = NULL
fpr.wave = NULL
uni.p.wave = unique(p.wave)
for(i in 1:length(uni.p.wave)){
    #i = 1
    wh = which(p.wave <= uni.p.wave[i])
    sig.wave[i] = length(wh)
    fdp.wave[i] = 1 - (sum(d.wave[wh])/length(wh))
    tpr.wave[i] = sum(d.wave[wh])/574
    fpr.wave[i] = (length(wh) - sum(d.wave[wh]))/574 
}



given.FDR = 0.05
posi = max(which(fdp.wave <= given.FDR))
tpr.wave[posi]
# 0.93
fdp.wave[posi]
# 0.049


## for ms
pval.ms = as.numeric(c(pval.null.ms, pval.alt.ms))
disc.ms = c(rep(0,574), rep(1,574))
rnk.ms = order(pval.ms)
p.ms = pval.ms[rnk.ms]
d.ms = disc.ms[rnk.ms]


fdp.ms = NULL
sig.ms = NULL
tpr.ms = NULL
fpr.ms = NULL
uni.p.ms = unique(p.ms)
for(i in 1:length(uni.p.ms)){
    wh = which(p.ms <= uni.p.ms[i])
    sig.ms[i] = length(wh)
    fdp.ms[i] = 1 - (sum(d.ms[wh])/length(wh))
    tpr.ms[i] = sum(d.ms[wh])/574
    fpr.ms[i] = (length(wh) - sum(d.ms[wh]))/574

}

given.FDR = 0.05
posi = max(which(fdp.ms <= given.FDR))
tpr.ms[posi]
# 0.86
fdp.ms[posi]
# 0.049




## for JAR
pval.JAR = as.numeric(c(pval.null.JAR, pval.alt.JAR))
disc.JAR = c(rep(0,574), rep(1,574))
rnk.JAR = order(pval.JAR)
p.JAR = pval.JAR[rnk.JAR]
d.JAR = disc.JAR[rnk.JAR]


fdp.JAR = NULL
sig.JAR = NULL
tpr.JAR = NULL
fpr.JAR = NULL
uni.p.JAR = unique(p.JAR)
for(i in 1:length(uni.p.JAR)){
    wh = which(p.JAR <= uni.p.JAR[i])
    sig.JAR[i] = length(wh)
    fdp.JAR[i] = 1 - (sum(d.JAR[wh])/length(wh))
    tpr.JAR[i] = sum(d.JAR[wh])/574
    fpr.JAR[i] = (length(wh) - sum(d.JAR[wh]))/574

}

given.FDR = 0.05
posi = max(which(fdp.JAR <= given.FDR))
tpr.JAR[posi]
# 0.8536
fdp.JAR[posi]
# 0.048





pdf("simu_574_ROC_wave_ms_JAR_80.pdf")

xmax = 1
ymax = 1
plot(c(0, fpr.wave), c(0, tpr.wave), xlim = c(0, xmax), ylim = c(0, ymax), xlab = "FPR (1 - Specificity)", ylab ="TPR (Sensitivity)", col = "red", type="l") 
points(c(0, fpr.ms), c(0, tpr.ms), type="l", col ="blue")
points(c(0, fpr.JAR), c(0, tpr.JAR), type="l", col ="green")

legend(0.54, 0.2, c("Wavelet based approach", "100bp window approach", "multi-scale") , col = c("red", "green", "blue"), lty = c(1,1), text.col = "black",merge = FALSE, bg = "white")
dev.off()









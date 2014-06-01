

##### read pval 1000 ######

##### Null ms ######
load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/null/multiscale/sum/pval.Robj")
pval.null.ms = as.numeric(pval_list)
done.null.ms = done_res


sum(done.null.ms)
# 100

max(pval.null.ms, na.rm = TRUE)
min(pval.null.ms, na.rm = TRUE)


##### alt ms  ######
load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/multiscale/sum/pval.Robj")
pval.alt.ms = as.numeric(pval_list)
done.alt.ms = done_res


sum(done.alt.ms)
# 100

max(pval.alt.ms, na.rm = TRUE)
min(pval.alt.ms, na.rm = TRUE)


delIX = c(305, 368, 451, 464)


pval.null.ms = pval.null.ms[-delIX]
pval.alt.ms = pval.alt.ms[-delIX]


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
    #i = 1
    wh = which(p.ms <= uni.p.ms[i])
    sig.ms[i] = length(wh)
    fdp.ms[i] = 1 - (sum(d.ms[wh])/length(wh))
    tpr.ms[i] = sum(d.ms[wh])/574
    fpr.ms[i] = (length(wh) - sum(d.ms[wh]))/574 
}


#tpr.ms[88:89]
given.FDR = 0.05
posi = max(which(fdp.ms <= given.FDR))
tpr.ms[posi]
# 0.98
fdp.ms[posi]
# 0.0485

fpr.ms.1000 = fpr.ms
tpr.ms.1000 = tpr.ms

pval.ms = as.numeric(c(pval.null.ms, pval.alt.ms))
disc.ms = c(rep(0,574), rep(1,574))
rank.ms = rank(pval.ms)
rank.ms.1000 = rank.ms




##### read pval 10000 ######

##### Null ms ######
load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/null/multiscale/sum_10000/pval.Robj")
pval.null.ms = as.numeric(pval_list)
done.null.ms = done_res


sum(done.null.ms)
# 100

max(pval.null.ms, na.rm = TRUE)
min(pval.null.ms, na.rm = TRUE)


##### alt ms  ######
load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/multiscale/sum_10000/pval.Robj")
pval.alt.ms = as.numeric(pval_list)
done.alt.ms = done_res


sum(done.alt.ms)
# 100

max(pval.alt.ms, na.rm = TRUE)
min(pval.alt.ms, na.rm = TRUE)


delIX = c(305, 368, 451, 464)


pval.null.ms = pval.null.ms[-delIX]
pval.alt.ms = pval.alt.ms[-delIX]


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
    #i = 1
    wh = which(p.ms <= uni.p.ms[i])
    sig.ms[i] = length(wh)
    fdp.ms[i] = 1 - (sum(d.ms[wh])/length(wh))
    tpr.ms[i] = sum(d.ms[wh])/574
    fpr.ms[i] = (length(wh) - sum(d.ms[wh]))/574 
}


#tpr.ms[88:89]
given.FDR = 0.05
posi = max(which(fdp.ms <= given.FDR))
tpr.ms[posi]
# 0.98
fdp.ms[posi]
# 0.0485

fpr.ms.10000 = fpr.ms
tpr.ms.10000 = tpr.ms

pval.ms = as.numeric(c(pval.null.ms, pval.alt.ms))
disc.ms = c(rep(0,574), rep(1,574))
rank.ms = rank(pval.ms)
rank.ms.10000 = rank.ms




pdf("simu574_ROCandRank_80_21_1000_vs_10000_ms.pdf")

xmax = 1
ymax = 1
plot(c(0, fpr.ms.10000), c(0, tpr.ms.10000), xlim = c(0, xmax), ylim = c(0, ymax), xlab = "FPR (1 - Specificity)", ylab ="TPR (Sensitivity)", col = "red", type="l") 
points(c(0, fpr.ms.1000), c(0, tpr.ms.1000), type="l", col ="blue")
legend(0.4, 0.12, c("multi-scale 10000 permutatons", "multi-scale 1000 permutation") , col = c("red", "blue"), lty = c(1,1), text.col = "black",merge = FALSE, bg = "white")




status.info = disc.ms

rank.10000 = rank.ms.10000
rank.1000 = rank.ms.1000

ymax = max(rank.1000)
ymin = min(rank.1000)
xmax = max(rank.10000)
plot(rank.10000, rank.1000, xlim = c(0, xmax), ylim = c(0, ymax), xlab = "rank in 10000 permutations", ylab ="rank in 1000 permutations", col = c("green", "darkgreen")[status.info+1], type="p", main ="multi-scale based test (dark green: alt)")


dev.off()












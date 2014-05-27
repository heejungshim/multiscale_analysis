


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

#load("/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/multiscale/sum/pval.Robj")
#pval.alt.ms = as.numeric(pval_list)
#done.alt.ms = done_res



pdf("hist_null_578.pdf")
par(mfrow = c(3,1))
hist(pval.null.wave, main="null Wavelet", breaks = 47)
hist(pval.null.JAR, main="null 100-window", breaks = 47)
hist(pval.null.ms, main="null multi-scale", breaks = 47)

dev.off()


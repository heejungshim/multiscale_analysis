
## read information on 578 sites
data = read.table("/mnt/lustre/home/shim/wavelets/revision/etc/simu.578.sites.txt", header=T)
#names(data)
#[1] "index"       "chr"         "site"        "genoIX"      "FDR.10.wave"
wh = which(is.na(data$genoIX) == TRUE)
data[wh,]

#    index chr site genoIX FDR.10.wave
#305 28769  12   42     NA           1
#368 34733  15  856     NA           1
#464 43187  19  750     NA           1

# let's get genoIX for those dsQTLs.
logLRs = read.table("/mnt/lustre/home/shim/wavelets/data/DNase/run/wavelets_01_step1_mean2_new_backup/output/SIM.12.42.fph.BFs.txt", as.is = TRUE)
dim(logLRs)
# 19 1026
which.max(as.numeric(logLRs[,2]))
# 17

logLRs = read.table("/mnt/lustre/home/shim/wavelets/data/DNase/run/wavelets_01_step1_mean2_new_backup/output/SIM.15.856.fph.BFs.txt", as.is = TRUE)
dim(logLRs)
# 12 1026
which.max(as.numeric(logLRs[,2]))
# 6


logLRs = read.table("/mnt/lustre/home/shim/wavelets/data/DNase/run/wavelets_01_step1_mean2_new_backup/output/SIM.19.750.fph.BFs.txt", as.is = TRUE)
dim(logLRs)
# 28 1026
which.max(as.numeric(logLRs[,2]))
# 3


genoIX.temp = c(17,6,3)
data$genoIX[wh] = genoIX.temp
sum(is.na(data$genoIX))

write.table(data, file = "/mnt/lustre/home/shim/wavelets/revision/etc/simu.578.sites.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)


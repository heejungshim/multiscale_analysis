###############################################
# implement function to remove PCR artifacts
###############################################

# input matrix
set.seed(1)
data = matrix(data=rpois(150*4, 0.1), nc = 150, nr = 4)
data[1,70] = 500
data[3,110] = 500
data[2,30] = 500
data[3,30] = 500


max.val = apply(data, 2, max)
candidate.posi = which(max.val > 1)
candidate.posi
#[1] 30 70  81 110 137
max.val[candidate.posi]
#[1] 500 500   2 500   2


# inside of function
#win.half.size = 50
#prop.thresh = 0.9
max.val = apply(data, 2, max)
candidate.posi = which(max.val > 1)


get.pcr.artifacts.posi <- function(posi, data, win.half.size = 50, prop.thresh = 0.9){
   
    st.win = max(1, posi - win.half.size)        
    en.win = st.win + win.half.size*2
    if(en.win > dim(data)[2]){
        en.win = dim(data)[2]
        st.win = en.win - win.half.size*2
    }
    prop = data[,posi]/apply(data[,st.win:en.win], 1, sum)
    pcr.posi = which(prop > prop.thresh)
    return(pcr.posi)
}



num.sam = dim(data)[1]
len = length(candidate.posi)
if(len > 0){
    pcr.posi.list = lapply(candidate.posi, get.pcr.artifacts.posi, data = data, win.half.size = 50, prop.thresh = 0.9)
    pcr.posi = which(sapply(pcr.posi.list, length) > 0)
    len.pcr.posi = length(pcr.posi)
    if(len.pcr.posi > 0){
        for(p in 1:len.pcr.posi){
            pcr.sam = pcr.posi.list[[pcr.posi[p]]]
            ix = candidate.posi[pcr.posi[p]]
            if(length(pcr.sam) == num.sam){
                data[,ix] = 1
            }else{
                data[, ix] = max(1, ceiling(mean(data[-pcr.sam,ix])))
            }
                
}


max.val = apply(data, 2, max)
candidate.posi = which(max.val > 1)
candidate.posi
#[1] 81 137
max.val[candidate.posi]
#[1] 2 2




## check implemention in my.utils.R


setwd("/mnt/lustre/home/shim/multiscale_analysis")

library("multiseq")
library("ashr")


multiscale.analysis.repodir <- scan(".multiscale_analysis.repodir.txt", what=character())
source(paste0(multiscale.analysis.repodir, "/src/R/my.utils.R"))
        
## input matrix
set.seed(1)
data = matrix(data=rpois(150*4, 0.1), nc = 150, nr = 4)
data[1,70] = 500
data[3,110] = 500
data[2,30] = 500
data[3,30] = 500


max.val = apply(data, 2, max)
candidate.posi = which(max.val > 1)
candidate.posi
#[1] 30 70  81 110 137
max.val[candidate.posi]
#[1] 500 500   2 500   2
        
source(paste0(multiscale.analysis.repodir, "/src/R/my.utils.R"))
res = remove.pcr.artifacts(data)

res$posi.with.pcr.artifacts
# 30 70 110 
data = res$data
max.val = apply(data, 2, max)
candidate.posi = which(max.val > 1)
candidate.posi
#[1] 81 137
max.val[candidate.posi]
#[1] 2 2




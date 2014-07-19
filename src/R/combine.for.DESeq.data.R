## `combine.for.DESeq.data.R' contains scrits to read data (prepread by `preprocess.wave.DESeq.roger.ATACseq.R'), make them as one matrix (as DESeq input format), and save them as a file. 
##
## Example Usage (see command in /mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/com) : R CMD BATCH --no-save --no-restore "--args chr=$SGE_TASK_ID st.sites=1 en.sites=10 wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/' siteSize=2048 treatment='Copper' null=FALSE strand='plus' window.size=100 numSam=6" /mnt/lustre/home/shim/multiscale_analysis/src/R/combine.for.DESeq.data.R
##
##
## chr : chromosome
## st.sites : first site to combine
## en.sites :  last site to combine 
## wd.path : working directory path
## siteSize : site size
## treatment : treatment name
## null : indicate whether it's null (control 1 vs control 2) or alternative data
## strand : 'both', 'plus', 'minus'; add two strands, use + strand, or use - strand
## window.size : window size we consider for DESeq analysis
## numSam : number of samples 
##
##
## Copyright (C) 2014 Heejung Shim
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.



#chr=1
#st.sites=1
#en.sites=10
#wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/deseq/'
#siteSize=2048
#treatment='Copper'
#null=FALSE
#strand='both'
#strand='plus'
#strand='minus'
#window.size=100
#numSam = 6


args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))
eval(parse(text=args[[2]]))
eval(parse(text=args[[3]]))
eval(parse(text=args[[4]]))
eval(parse(text=args[[5]]))
eval(parse(text=args[[6]]))
eval(parse(text=args[[7]]))
eval(parse(text=args[[8]]))
eval(parse(text=args[[9]]))
eval(parse(text=args[[10]]))




## set up working directory 
setwd(wd.path)

## make directory name  
if(!null){
    dir.name=paste0(treatment,".", siteSize, ".", strand,  ".", window.size, ".alt")
}else{
    dir.name=paste0(treatment,".", siteSize, ".", strand, ".", window.size, ".null")
}

## make input/output directory name and make output directory 
input.dir.path = paste0(wd.path, dir.name, ".data/") 
output.dir.path = paste0(wd.path, dir.name, ".run/") 
if(!file.exists(output.dir.path)){
    dir.create(output.dir.path)
}


numSites = en.sites - st.sites + 1
numC = siteSize%/%window.size
numRow = numSites*numC

st.list = (1:numSites)*numC - (numC -1)
en.list = (1:numSites)*numC

DESeq.dat = matrix(data=NA, nc = numSam, nr = numRow)
done = rep(NA, numSites)

for(sites in 1:numSites){

    this.path = paste0(input.dir.path, "data.", chr, ".", sites, ".txt")
    if(file.exists(this.path)== FALSE){		
        done[sites] = FALSE
    }else{
        if(file.info(this.path)$size == 0){
            done[sites] = FALSE
        }else{
            dat = read.table(this.path)
            if(sum(dim(dat)==c(numC, numSam))==2){
                DESeq.dat[st.list[sites]:en.list[sites],] = as.matrix(dat)
                done[sites] = TRUE
            }else{
                done[sites] = FALSE
            }
        }
    }
}


this.out.path = paste0(output.dir.path, "data.", chr, ".txt")
write.table(DESeq.dat, file=this.out.path, quote=FALSE, row.names = FALSE, col.names = FALSE)


this.out.path = paste0(output.dir.path, "err.", chr, ".txt")
write.table(which(done==FALSE), file=this.out.path, quote=FALSE, row.names = FALSE, col.names = FALSE)

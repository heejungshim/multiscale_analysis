## `run.DESeq2.on.simulation.R' contains scrits to run DESeq2 on simulated data.
## 
##
## Example Usage (see command in /mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint_new/sum/DESeq/com/): R CMD BATCH --no-save --no-restore "--args wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint_new/' numSam=6 read.depth.ratio=1 num.simu=500 siteSize=1024 filter.cut=0" /mnt/lustre/home/shim/multiscale_analysis/src/R/run.DESeq2.on.simulation.R
##
##
## wd.path : working directory path
## numSam : number of sample
## read.depth.ratio : library read depth (either x0.5, x1, x2, x4)
## num.simu : number of simulations
## siteSize : site size in simulation 
## filter.cut : analysis includes window with read count > filter.cut
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


args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))
eval(parse(text=args[[2]]))
eval(parse(text=args[[3]]))
eval(parse(text=args[[4]]))
eval(parse(text=args[[5]]))
eval(parse(text=args[[6]]))


##wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint_new/'
##numSam = 6
##read.depth.ratio = 1
##num.simu = 500
##siteSize = 1024
##filter.cut = 0

library("DESeq2")

window.size.list = c(100, 300, 1024)

## set up working directory 
setwd(wd.path)

## set up directory name
if(read.depth.ratio==0.5){
  dir.name=paste0("halfread.", numSam, "ind")
}
if(read.depth.ratio==1){
  dir.name=paste0("fullread.", numSam, "ind")
}
if(read.depth.ratio==2){
  dir.name=paste0("fullread2.", numSam, "ind")
}
if(read.depth.ratio==4){
  dir.name=paste0("fullread4.", numSam, "ind")
}

## set up name list to run DESeq
genoD = c(rep("T", numSam/2), rep("C", numSam/2))
name.list = genoD
for(i in 1:floor(numSam/2)){
  name.list[i] = paste0(genoD[i], i)
  name.list[numSam/2 + i] = paste0(genoD[numSam/2 + i], i)
}

for(ww in 1:3){

  window.size = window.size.list[ww]
  numC = siteSize%/%window.size

  ## read null data
  input.path = paste0(wd.path, "null/DESeq/", dir.name, ".", window.size, ".run/")
  deseq.data.null = read.table(paste0(input.path, "data.txt"))

  ## read alt data
  input.path = paste0(wd.path, "alt/DESeq/", dir.name, ".", window.size, ".run/")
  deseq.data.alt = read.table(paste0(input.path, "data.txt"))

  ## combine data
  deseq.data = rbind(deseq.data.null, deseq.data.alt)

  ## Converts the colData table to a dataframe with vector labels
  colData <- data.frame(row.names = as.vector(name.list), t = as.vector(genoD))

  ## filter data
  rsum = rowSums (deseq.data)
  use = ((rsum > filter.cut) & (!is.na(rsum)))
  countData.filtered = deseq.data[use, ]

  ## Perform DESeq2 Analysis
  ddsTvC <- DESeqDataSetFromMatrix(
      countData=countData.filtered,
      colData=colData,
      design = ~ t)
  dds <- DESeq(ddsTvC)
  res <- results(dds)

  ## get p-value
  pval.vec = rep(NA, length(use))
  pval.vec[use==TRUE] = res$pvalue
  pval.filtered =matrix(pval.vec,ncol=numC,byrow=T)

  ## get minimum p-value for each site
  min.pval=apply(pval.filtered,1,min,na.rm=TRUE)
  min.pval[is.infinite(min.pval)] = NA

  ## try to save output
  output.path = paste0(wd.path, "sum/DESeq/", dir.name, ".", window.size)

  ## output min.pval
  write.table(min.pval, file = paste0(output.path, ".min.pval.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)

  ## output object just in case we want to take a close look!
  save("res", "use", file =paste0(output.path, ".Robj"))

}

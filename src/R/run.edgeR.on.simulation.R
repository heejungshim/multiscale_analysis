## `run.edgeR.on.simulation.R' contains scrits to run edgeR on simulated data.
## 
##
## Example Usage (see command in /mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_manyQTLfinal_v1/sum/edgeR/com/): R CMD BATCH --no-save --no-restore "--args wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_manyQTLfinal_v1/' numSam=10 read.depth.ratio=1 num.simu=578 siteSize=1024 filter.cut=0 combineNullandAlt=TRUE separateNullandAlt=TRUE" /mnt/lustre/home/shim/multiscale_analysis/src/R/run.edgeR.on.simulation.R
##
##
## wd.path : working directory path
## numSam : number of sample
## read.depth.ratio : library read depth (either x0.5, x1, x2, x4)
## num.simu : number of simulations
## siteSize : site size in simulation 
## filter.cut : analysis includes window with read count > filter.cut
## null.path : path to directory where null data have been saved.
## combineNullandAlt : if it is true, the script runs software after combining null and alternative data sets
## separateNullandAlt : if it is true, the script runs software after combining null and alternative data sets
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
eval(parse(text=args))

if(!exists("null.path")){
  null.path = NULL
}
if(!exists("combineNullandAlt")){
  combineNullandAlt = TRUE
}
if(!exists("separateNullandAlt")){
  separateNullandAlt = TRUE
}


##wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_manyQTLfinal_v2/'
##numSam=10
##read.depth.ratio=1
##num.simu=578
##siteSize=1024
##filter.cut=0
##null.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_manyQTLfinal_v1/'
##combineNullandAlt = TRUE
##separateNullandAlt = TRUE

library("edgeR")

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

## set up gp to run edgeR
gp = factor(c(rep(0,numSam/2), rep(1,numSam/2)) ,labels=c("A","B"))

for(ww in 1:3){

  window.size = window.size.list[ww]
  numC = siteSize%/%window.size

  ## read null data
  if(is.null(null.path)){
    input.path = paste0(wd.path, "null/DESeq/", dir.name, ".", window.size, ".run/")
  }else{
    input.path = paste0(null.path, "null/DESeq/", dir.name, ".", window.size, ".run/")
  }    
  edgeR.data.null = read.table(paste0(input.path, "data.txt"))

  ## read alt data
  input.path = paste0(wd.path, "alt/DESeq/", dir.name, ".", window.size, ".run/")
  edgeR.data.alt = read.table(paste0(input.path, "data.txt"))

  if(combineNullandAlt){
    ## combine data
    edgeR.data = rbind(edgeR.data.null, edgeR.data.alt)

    ## filter data
    rsum = rowSums(edgeR.data)
    use = ((rsum > filter.cut) & (!is.na(rsum)))
    countData.filtered = edgeR.data[ use, ]

    ## Perform edgeR Analysis
    res = DGEList(counts=countData.filtered, group=gp)
    res = calcNormFactors(res, method="RLE")
    res = estimateCommonDisp(res)
    res = estimateTagwiseDisp(res)
    res = exactTest(res, dispersion="auto")

    ## get p-value
    pval.vec = rep(NA, length(use))
    pval.vec[use==TRUE] = res$table$PValue
    pval.filtered =matrix(pval.vec,ncol=numC,byrow=T)

    ## get minimum p-value for each site
    min.pval=apply(pval.filtered,1,min,na.rm=TRUE)
    min.pval[is.infinite(min.pval)] = NA

    ## try to save output
    output.path = paste0(wd.path, "sum/edgeR/", dir.name, ".", window.size)

    ## output min.pval
    write.table(min.pval, file = paste0(output.path, ".min.pval.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)

    ## output object just in case we want to take a close look!
    save("res", "use", file =paste0(output.path, ".Robj"))

  }

  if(separateNullandAlt){

    ## for null data
    edgeR.data = edgeR.data.null

    ## filter data
    rsum = rowSums(edgeR.data)
    use = ((rsum > filter.cut) & (!is.na(rsum)))
    countData.filtered = edgeR.data[ use, ]

    ## Perform edgeR Analysis
    res = DGEList(counts=countData.filtered, group=gp)
    res = calcNormFactors(res, method="RLE")
    res = estimateCommonDisp(res)
    res = estimateTagwiseDisp(res)
    res = exactTest(res, dispersion="auto")

    ## get p-value
    pval.vec = rep(NA, length(use))
    pval.vec[use==TRUE] = res$table$PValue
    pval.filtered =matrix(pval.vec,ncol=numC,byrow=T)

    ## get minimum p-value for each site
    min.pval=apply(pval.filtered,1,min,na.rm=TRUE)
    min.pval[is.infinite(min.pval)] = NA

    min.pval.null = min.pval
    res.null = res
    use.null = use

    ## for alt data
    edgeR.data = edgeR.data.alt

    ## filter data
    rsum = rowSums(edgeR.data)
    use = ((rsum > filter.cut) & (!is.na(rsum)))
    countData.filtered = edgeR.data[ use, ]

    ## Perform edgeR Analysis
    res = DGEList(counts=countData.filtered, group=gp)
    res = calcNormFactors(res, method="RLE")
    res = estimateCommonDisp(res)
    res = estimateTagwiseDisp(res)
    res = exactTest(res, dispersion="auto")

    ## get p-value
    pval.vec = rep(NA, length(use))
    pval.vec[use==TRUE] = res$table$PValue
    pval.filtered =matrix(pval.vec,ncol=numC,byrow=T)

    ## get minimum p-value for each site
    min.pval=apply(pval.filtered,1,min,na.rm=TRUE)
    min.pval[is.infinite(min.pval)] = NA

    min.pval.alt = min.pval
    res.alt = res
    use.alt = use

    min.pval = c(min.pval.null, min.pval.alt)
    use = c(use.null, use.alt)

    ## try to save output
    output.path = paste0(wd.path, "sum/edgeR/sep.", dir.name, ".", window.size)

    ## output min.pval
    write.table(min.pval, file = paste0(output.path, ".min.pval.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)

    ## output object just in case we want to take a close look!
    save("res.alt", "res.null", "use", "use.alt", "use.null", file =paste0(output.path, ".Robj"))
    
  }
  
}

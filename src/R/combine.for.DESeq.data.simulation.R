## `combine.for.DESeq.data.simulation.R' contains scrits to read simulated data, make them as one matrix (as DESeq2 input format), and save them as a file. 
##
## Example Usage (see command in /mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint_new/null/DESeq/com/) R CMD BATCH --no-save --no-restore "--args wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint_new/null/DESeq/' numSam=6 read.depth.ratio=1 num.simu=500 siteSize=1024" /mnt/lustre/home/shim/multiscale_analysis/src/R/combine.for.DESeq.data.simulation.R
##
##
## wd.path : working directory path
## numSam : number of samples
## read.depth.ratio : library read depth ratio (either 0.5, 1, 2, 4)
## num.simu : number of simulations
## siteSize : site size in simulation 
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



##wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint_new/null/DESeq/'
##numSam = 6
##read.depth.ratio = 1
##num.simu = 500
##siteSize = 1024


args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))
eval(parse(text=args[[2]]))
eval(parse(text=args[[3]]))
eval(parse(text=args[[4]]))
eval(parse(text=args[[5]]))


window.size.list = c(100, 300, 1024)

## set up working directory 
setwd(wd.path)

## set up directory name
if(read.depth.ratio < 2){
  if(read.depth.ratio==0.5){
    dir.name=paste0("halfread.", numSam, "ind")
  }
  if(read.depth.ratio==1){
    dir.name=paste0("fullread.", numSam, "ind")
  }
}else{
  dir.name=paste0("fullread", read.depth.ratio, ".", numSam, "ind")
}

for(ww in 1:3){

  window.size = window.size.list[ww]
  
  ## make input/output directory name and make output directory 
  input.dir.path = paste0(wd.path, dir.name, ".", window.size, ".data/") 
  output.dir.path = paste0(wd.path, dir.name, ".", window.size, ".run/") 
  if(!file.exists(output.dir.path)){
    dir.create(output.dir.path)
  }

  numC = siteSize%/%window.size
  numRow = num.simu*numC
  
  st.list = (1:num.simu)*numC - (numC -1)
  en.list = (1:num.simu)*numC
  
  DESeq.dat = matrix(data=NA, nc = numSam, nr = numRow)
  done = rep(NA, num.simu)
  
  for(sites in 1:num.simu){
    
    this.path = paste0(input.dir.path, "DESeq.", sites, ".txt")
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
  
  this.out.path = paste0(output.dir.path, "data.txt")
  write.table(DESeq.dat, file=this.out.path, quote=FALSE, row.names = FALSE, col.names = FALSE)

  this.out.path = paste0(output.dir.path, "err.txt")
  write.table(which(done==FALSE), file=this.out.path, quote=FALSE, row.names = FALSE, col.names = FALSE)

}


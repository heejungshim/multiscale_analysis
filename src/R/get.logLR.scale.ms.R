## `get.logLR.scale.ms.R' contains scrits to collect logLR for overall mean and logLR for shape from multiseq analysis, save them as a vector, and output them as R object.
##
## Example Usage (see command in /mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/multiscale/com/Copper.2048.plus.alt.post) : /data/tools/R-3.1.0/bin/R CMD BATCH --no-save --no-restore "--args chr=10 st.sites=1 en.sites=7205 path.output.dir='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/multiscale/Copper.2048.plus.alt.output/' path.sum.dir='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/multiscale/Copper.2048.plus.alt.sum/' output.file.name='Copper.2048.plus.alt'" /mnt/lustre/home/shim/multiscale_analysis/src/R/get.logLR.scale.ms.R
##
##
## chr : chromosome (if chr=NULL, we assume there is no chromosome). Otherwise, output file format is chr.sites
## st.sites : first site to combine
## en.sites :  last site to combine 
## path.output.dir : path to directory which contain logLR files
## path.sum.dir : path to directory where we will save R object
## output.file.name : R object name
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
#en.sites=100
#path.output.dir='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/multiscale/Copper.2048.plus.alt.output/'
#path.sum.dir='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/run/multiscale/Copper.2048.plus.alt.sum/'
#output.file.name='Copper.2048.plus.alt'

args = (commandArgs(TRUE))
eval(parse(text=args))


numSites = en.sites - st.sites + 1
done = rep(NA, numSites)
logLR.mean = rep(NA, numSites)
logLR = rep(NA, numSites)


path = paste0(path.output.dir, "res.")

for(i in st.sites:en.sites){

    if(is.null(chr)){
        path.each  = paste0(path, i, ".out")
    }else{
        path.each  = paste0(path, chr, ".", i, ".out")
    }

    if(file.exists(path.each)== FALSE){		
        done[i] = FALSE
    }else{
        if(file.info(path.each)$size == 0){
            done[i] = FALSE
        }else{			
            dat = scan(path.each, what=double(), quiet=TRUE)
            llr.mean = dat[length(dat)]
            llr = sum(dat[2:(length(dat)-1)])
            if((!is.na(llr.mean)) & (!is.na(llr))){
                logLR.mean[i] = llr.mean
                logLR[i] = llr   
                done[i] = TRUE
            }else{
                done[i] = FALSE
            }
        }
    }
}

## make summary directory
if(!file.exists(path.sum.dir)){
    dir.create(path.sum.dir)
}

if(is.null(chr)){
    save("logLR", "logLR.mean", "done", file =paste0(path.sum.dir, output.file.name, ".scale.Robj"))
}else{
    save("logLR", "logLR.mean", "done", file =paste0(path.sum.dir, output.file.name, ".", chr, ".scale.Robj"))
}







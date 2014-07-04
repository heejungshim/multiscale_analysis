#!/usr/bin/env Rscript

## Aim : This file contains Rscripts to get effect sizes and raw data from 578 dsQTLs idetinfied by Shim and Stephens 2014 (see its Supplementary Materials for details of those dsQTLs. This information will be used in simulation.
## I modified two scripts: "/mnt/lustre/home/shim/wavelets/revision/code/simulation.explore.578.sites.R" and "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/code/better.effect.size.R"
##
## Usage R CMD BATCH --no-save --no-restore "--args ss=$SGE_TASK_ID" prepare.data.for.simulation.R
## See /mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/com/com.effect.size.sh
##
## 
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+

## ss = 305
args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))


setwd("/mnt/lustre/home/shim/multiscale_analysis")
multiscale.analysis.repodir <- scan(".multiscale_analysis.repodir.txt", what=character())
WaveQTL.repodir <- scan(".WaveQTL.repodir.txt", what=character())

## set working directory
wd.path = paste0(multiscale.analysis.repodir, "/analysis/simulation/sample_size/simulation_578/code/")
setwd(wd.path)


## set path to files 

## Path to directory which contain DNase-seq data as hdf5 format, 
hdf5.data.path = "/mnt/lustre/data/internal/genome_db/hg18/dnase/"
## Path to mappability information as hdf5 format
hdf5.mapp.path  = "/mnt/lustre/data/internal/genome_db/hg18/mappability/roger_20bp_mapping_uniqueness.h5"

## path to directory which contains information on SNPs located region of interest
geno.info.dir.path = "/mnt/lustre/home/shim/wavelets/data/DNase/geno_01_step1/geno/"

## path to directory which contains genotype data
geno.dir.path = "/mnt/lustre/home/shim/wavelets/data/DNase/geno_01_step1/geno_maf/"


## path to directory which contains location information on 578 sites
locus.path = "/mnt/lustre/home/shim/wavelets/data/DNase/region_01_sel_step1/"

## path to output directory
output.path = paste0(multiscale.analysis.repodir, "/analysis/simulation/sample_size/simulation_578/data/")

## read a list of individual IDs
inds.IDs = scan(paste0(WaveQTL.repodir, "/data/Shim_2014_etc/DNaseI.individuals.oneline.txt"), what="")

## read functions to read DNase data and preprocess 
source(paste0(multiscale.analysis.repodir, "/src/R/prepare.DNase.funcs.R"))

## read information on 578 sites
data = read.table("/mnt/lustre/home/shim/wavelets/revision/etc/simu.578.sites.txt", header=T)
#names(data)
#[1] "index"       "chr"         "site"        "genoIX"      "FDR.10.wave"


chr.list = data$chr
site.list = data$site
genoIX.list = data$genoIX


## for each dsQTL, let's prepare data

#for(ss in 1:578){


    chr = chr.list[ss]
    site = site.list[ss]
    genoIX = genoIX.list[ss]

    ## read location information 
    path = paste0(locus.path, "chr", chr, ".loc")
    loc_dat = read.table(path, as.is = TRUE)
    chrIX = loc_dat[site,1]
    locus.start = loc_dat[site,2]
    locus.end = loc_dat[site, 3] - 1

    ## path to genotype information (to correct for cutting preference)
    geno.info.path = paste0(geno.info.dir.path, "maf_chr", chr, ".", site, ".geno")

    ## run function to read DNase data 
    res = read.DNase.data(hdf5.data.path = hdf5.data.path, hdf5.mapp.path = hdf5.mapp.path, geno.info.path = geno.info.path, inds.IDs = inds.IDs, chrIX = chrIX, locus.start = locus.start  , locus.end = locus.end)

    phenoD = res$DNase.dat

    ## get genotype informaiton
    geno.path = paste0(geno.dir.path, "chr", chr, ".", site, ".geno")
    genoF = read.table(geno.path, as.is = TRUE)
    genoD = genoF[genoIX, 4:73]
    genoR = as.numeric(round(genoD))

    ## Estimate effect size 
    wh0 = which(genoR == 0)
    wh1 = which(genoR == 1)

    if(min(length(wh0), length(wh1)) == 0){
        path.output = paste0(output.path, "warnings.", ss)
        cat("no individauls either one of groups", file = path.output)
    }else{
       
        ## Take a mean profile
        sig0 = apply(phenoD[wh0,], 2, mean)
        sig1 = apply(phenoD[wh1,], 2, mean)

        ## handle 0 count
        wh.zero = which(sig0 == 0)
        if(length(wh.zero) > 0){ 
            sig0[wh.zero] = 1/70
        }
        wh.zero = which(sig1 == 0)
        if(length(wh.zero) > 0){ 
            sig1[wh.zero] = 1/70
        }

        ## denoise mean curves using wavethresh
        library(wavethresh)

        w.sig0 = wd(sig0)
        w.sig0.smooth = threshold(w.sig0)
        sig0.smooth = wr(w.sig0.smooth)

        w.sig1 = wd(sig1)
        w.sig1.smooth = threshold(w.sig1)
        sig1.smooth = wr(w.sig1.smooth)

        ## get smooth signal as a ratio
        smooth.ratio = sig1.smooth/sig0.smooth

        ## output information to files
        path.output = paste0(output.path, "raw.dat.", ss)
        write.table(phenoD, path.output, row.names = FALSE, col.names = FALSE, quote = FALSE)

        path.output = paste0(output.path, "smooth.ratio.", ss)
        write.table(smooth.ratio, path.output, row.names = FALSE, col.names = FALSE, quote = FALSE)

    }
#}




    

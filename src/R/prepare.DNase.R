## `prepare.DNase.R' shows how to use functions in prepare.DNase.funcs.R
##
##
## Copyright (C) 2014 Heejung Shim
##
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



setwd("/mnt/lustre/home/shim/multiscale_analysis")


multiscale.analysis.repodir <- scan(".multiscale_analysis.repodir.txt", what=character())
WaveQTL.repodir <- scan(".WaveQTL.repodir.txt", what=character())
source(paste0(multiscale.analysis.repodir, "/src/R/prepare.DNase.funcs.R"))
source(paste0(multiscale.analysis.repodir, "/src/R/utils.R"))


#########################################
#### Users need to make changes in paths
#########################################

## Path to directory which contain DNase-seq data as hdf5 format, 
hdf5.data.path = "/mnt/lustre/data/internal/genome_db/hg18/dnase/"
## Path to mappability information as hdf5 format
hdf5.mapp.path  = "/mnt/lustre/data/internal/genome_db/hg18/mappability/roger_20bp_mapping_uniqueness.h5"

## list of individual IDs
inds.IDs = scan(paste0(WaveQTL.repodir, "/data/Shim_2014_etc/DNaseI.individuals.oneline.txt"), what="")

## location informaiton 
chrIX = "chr17"
locus.start = 10160989
locus.end = 10162013 - 1 

## genotype information (to correct for cutting preference)
geno.info.path = paste0(WaveQTL.repodir, "/data/dsQTL/chr17.10160989.10162012.2kb.cis.noMAFfilter.info")



res = read.DNase.data(hdf5.data.path = hdf5.data.path, hdf5.mapp.path = hdf5.mapp.path, geno.info.path = geno.info.path, inds.IDs = inds.IDs, chrIX = chrIX, locus.start = locus.start  , locus.end = locus.end)


str(res)
#List of 2
# $ DNase.dat  : num [1:70, 1:1024] 0 0 0 0 0 0 0 0 0 0 ...
# $ mappability: num [1:1024] 1 1 1 1 1 1 1 1 1 1 ...
dim(res$DNase.dat)
#[1]   70 1024
length(res$mappability)
#[1] 1024

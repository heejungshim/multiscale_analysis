#!/usr/bin/env Rscript

## Aim : This file contains Rscripts to select top 5% of hypersensitive 300bp sites (following a procedure Roger used) 
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+


# path to directory which contain information Roger provided 
folder = "/mnt/gluster/data/external_private_supp/roger_atacseq/deseq/"

# read genometile information
genometile <- read.table(paste(folder,"hg19.tiles.w300.s300.bed.gz",sep=""), sep="\t", as.is=T)
dim(genometile)
#[1] 10457248        3
genometile[1:2,]
#    V1  V2  V3
#1 chr1   0 300
#2 chr1 300 600


########################################
# For Copper vs Control 2
# 65 and 6500
########################################

tbarcode="N702"
cbarcode="N706"

rep1barcode="N501"
rep2barcode="N502"
rep3barcode="N503"

## Get read count for each tile
tRep1 <- read.table(paste(folder,tbarcode,rep1barcode,".bed.gz", sep=""), as.is=T)[,4]
tRep2 <- read.table(paste(folder,tbarcode,rep2barcode,".bed.gz", sep=""), as.is=T)[,4]
tRep3 <- read.table(paste(folder,tbarcode,rep3barcode,".bed.gz", sep=""), as.is=T)[,4]
cRep1 <- read.table(paste(folder,cbarcode,rep1barcode,".bed.gz", sep=""), as.is=T)[,4]
cRep2 <- read.table(paste(folder,cbarcode,rep2barcode,".bed.gz", sep=""), as.is=T)[,4]
cRep3 <- read.table(paste(folder,cbarcode,rep3barcode,".bed.gz", sep=""), as.is=T)[,4]

## Prepare Data Tables for DESeq
countData <- cbind(tRep1,tRep2,tRep3,cRep1,cRep2,cRep3)

## save total sum
rs.Copper <- rowSums(countData)



########################################
# For Selenium vs Control 2
# 65 and 6500
########################################

tbarcode="N703"
cbarcode="N706"

rep1barcode="N501"
rep2barcode="N502"
rep3barcode="N503"

## Get read count for each tile
tRep1 <- read.table(paste(folder,tbarcode,rep1barcode,".bed.gz", sep=""), as.is=T)[,4]
tRep2 <- read.table(paste(folder,tbarcode,rep2barcode,".bed.gz", sep=""), as.is=T)[,4]
tRep3 <- read.table(paste(folder,tbarcode,rep3barcode,".bed.gz", sep=""), as.is=T)[,4]
cRep1 <- read.table(paste(folder,cbarcode,rep1barcode,".bed.gz", sep=""), as.is=T)[,4]
cRep2 <- read.table(paste(folder,cbarcode,rep2barcode,".bed.gz", sep=""), as.is=T)[,4]
cRep3 <- read.table(paste(folder,cbarcode,rep3barcode,".bed.gz", sep=""), as.is=T)[,4]

## Prepare Data Tables for DESeq
countData <- cbind(tRep1,tRep2,tRep3,cRep1,cRep2,cRep3)

## save total sum
rs.Selenium <- rowSums(countData)




########################################
# For Retinoic vs Control 2
# 65 and 6500
########################################

tbarcode="N704"
cbarcode="N705"

rep1barcode="N501"
rep2barcode="N502"
rep3barcode="N503"

## Get read count for each tile
tRep1 <- read.table(paste(folder,tbarcode,rep1barcode,".bed.gz", sep=""), as.is=T)[,4]
tRep2 <- read.table(paste(folder,tbarcode,rep2barcode,".bed.gz", sep=""), as.is=T)[,4]
tRep3 <- read.table(paste(folder,tbarcode,rep3barcode,".bed.gz", sep=""), as.is=T)[,4]
cRep1 <- read.table(paste(folder,cbarcode,rep1barcode,".bed.gz", sep=""), as.is=T)[,4]
cRep2 <- read.table(paste(folder,cbarcode,rep2barcode,".bed.gz", sep=""), as.is=T)[,4]
cRep3 <- read.table(paste(folder,cbarcode,rep3barcode,".bed.gz", sep=""), as.is=T)[,4]

## Prepare Data Tables for DESeq
countData <- cbind(tRep1,tRep2,tRep3,cRep1,cRep2,cRep3)

## save total sum
rs.Retinoic <- rowSums(countData)



######################
# Select top 5%
######################
#rs.Retinoic
#rs.Selenium
#rs.Copper

rnk.Copper = rank(rs.Copper)
rnk.Selenium = rank(rs.Selenium)
rnk.Retinoic = rank(rs.Retinoic)

sum(rnk.Copper > 9950000)/length(rs.Copper)
#[1] 0.04740272
sum(rnk.Selenium > 9950000)/length(rs.Copper)
#[1] 0.04790591
sum(rnk.Retinoic > 9950000)/length(rs.Copper)
#[1] 0.0488531

min(rs.Copper[which(rnk.Copper > 9950000)])
#[1] 66
min(rs.Selenium[which(rnk.Selenium > 9950000)])
#[1] 61
min(rs.Retinoic[which(rnk.Retinoic > 9950000)])
#[1] 61


sum(rs.Copper >= 66)/length(rs.Copper)
#[1] 0.04740272
sum(rs.Copper >= 66)
#[1] 495702
sum(rs.Copper >= 6500)
#[1] 516
 
sum(rs.Selenium >= 61)/length(rs.Selenium)
#[1] 0.04790591
sum(rs.Selenium >= 61)
#[1] 500964
sum(rs.Selenium >= 6500)
#[1] 473
 
sum(rs.Retinoic >= 61)/length(rs.Retinoic)
#[1] 0.0488531
sum(rs.Retinoic >= 61)
#[1] 510869
sum(rs.Retinoic >= 6500)
#[1] 405


################################
# Save top 5% location information!
################################

sel.Copper = which(rs.Copper >= 66) 
sel.Selenium = which(rs.Selenium >= 61)
sel.Retinoic = which(rs.Retinoic >= 61)

#dim(genometile)
#total.len = length(rs.Copper)
#total.len

path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/info/"
write.table(genometile[sel.Copper,], file = paste0(path, "Copper.05.txt"), quote=FALSE,row.names = FALSE, col.names = FALSE)

write.table(genometile[sel.Selenium,], file = paste0(path, "Selenium.05.txt"), quote=FALSE,row.names = FALSE, col.names = FALSE)

write.table(genometile[sel.Retinoic,], file = paste0(path, "Retinoic.05.txt"), quote=FALSE,row.names = FALSE, col.names = FALSE)


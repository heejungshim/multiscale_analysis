#!/usr/bin/env Rscript

## Aim : This file contains Rscripts to select top 5% of hypersensitive 300bp sites (following a procedure Roger used) 
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+


# path to directory which contain autosome genometile information 
folder.process = "/mnt/gluster/data/external_private_supp/roger_atacseq/process_bams/"

# read genometile information
genometile <- read.table(paste(folder.process,"hg19.autosomes.tiles.w300.s300.bed.gz",sep=""), sep=" ", as.is=T)
dim(genometile)
#[1] 9603454       3
genometile[1:2,]
#    V1  V2  V3
#1 chr1   0 300
#2 chr1 300 600

index.autosome =scan(file = paste0(folder.process,"index.in.tiles.txt"))
length(index.autosome)
#[1] 9603454
index.autosome[1:2]
#[1] 1 2


# path to directory which contain information roger provided 
folder = "/mnt/gluster/data/external_private_supp/roger_atacseq/deseq/"


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
rs.Copper <- rowSums(countData)[index.autosome]


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
rs.Selenium <- rowSums(countData)[index.autosome]




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
rs.Retinoic <- rowSums(countData)[index.autosome]



######################
# Select top 5%
######################
#rs.Retinoic
#rs.Selenium
#rs.Copper

length(rs.Copper)*0.952
#[1] 9142488


rnk.Copper = rank(rs.Copper)
rnk.Selenium = rank(rs.Selenium)
rnk.Retinoic = rank(rs.Retinoic)



sum(rnk.Copper > 9142488)/length(rs.Copper)
#[1] 0.04743658
sum(rnk.Selenium > 9142488)/length(rs.Copper)
#[1] 0.04790641
sum(rnk.Retinoic > 9142488)/length(rs.Copper)
#[1] 0.04875329


min(rs.Copper[which(rnk.Copper > 9142488)])
#[1] 66
min(rs.Selenium[which(rnk.Selenium > 9142488)])
#[1] 61
min(rs.Retinoic[which(rnk.Retinoic > 9142488)])
#[1] 61

sum(rs.Copper >= 66)/length(rs.Copper)
#[1] 0.04743658
sum(rs.Copper >= 66)
#[1] 455555
  
sum(rs.Selenium >= 61)/length(rs.Selenium)
#[1] 0.04790641
sum(rs.Selenium >= 61)
#[1] 460067

sum(rs.Retinoic >= 61)/length(rs.Retinoic)
#[1] 0.04875329
sum(rs.Retinoic >= 61)
#[1] 468200
 


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


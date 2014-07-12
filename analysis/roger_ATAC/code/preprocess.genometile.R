#!/usr/bin/env Rscript

## Aim : This file contains Rscripts to take only autosome from genometile Roger provided and compute chromosome length
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


table(genometile[,1])
#                 chr1                 chr10                 chr11 
#               830836                451783                450022 
#chr11_gl000202_random                 chr12                 chr13 
#                  134                446173                383900 
#                chr14                 chr15                 chr16 
#               357832                341772                301183 
#                chr17       chr17_ctg5_hap1 chr17_gl000203_random 
#               270651                  5603                   125 
#chr17_gl000204_random chr17_gl000205_random chr17_gl000206_random 
#                  272                   582                   137 
#                chr18 chr18_gl000207_random                 chr19 
#               260258                    15                197097 
#chr19_gl000208_random chr19_gl000209_random  chr1_gl000191_random 
#                  309                   531                   355 
# chr1_gl000192_random                  chr2                 chr20 
#                 1825                810665                210086 
#                chr21 chr21_gl000210_random                 chr22 
#               160433                    93                171016 
#                 chr3                  chr4        chr4_ctg9_hap1 
#               660075                637181                  1969 
# chr4_gl000193_random  chr4_gl000194_random                  chr5 
#                  633                   639                603051 
#                 chr6         chr6_apd_hap1         chr6_cox_hap2 
#               570384                 15408                 15985 
#        chr6_dbb_hap3        chr6_mann_hap4         chr6_mcf_hap5 
#                15368                 15611                 16112 
#        chr6_qbl_hap6        chr6_ssto_hap7                  chr7 
#                15374                 16429                530463 
# chr7_gl000195_random                  chr8  chr8_gl000196_random 
#                  610                487881                   130 
# chr8_gl000197_random                  chr9  chr9_gl000198_random 
#                  124                470712                   301 
# chr9_gl000199_random  chr9_gl000200_random  chr9_gl000201_random 
#                  567                   624                   121 
#                 chrM        chrUn_gl000211        chrUn_gl000212 
#                   56                   556                   623 
#       chrUn_gl000213        chrUn_gl000214        chrUn_gl000215 
#                  548                   460                   576 
#       chrUn_gl000216        chrUn_gl000217        chrUn_gl000218 
#                  575                   574                   538 
#       chrUn_gl000219        chrUn_gl000220        chrUn_gl000221 
#                  598                   540                   518 
#       chrUn_gl000222        chrUn_gl000223        chrUn_gl000224 
#                  623                   602                   599 
#       chrUn_gl000225        chrUn_gl000226        chrUn_gl000227 
#                  704                    51                   428 
#       chrUn_gl000228        chrUn_gl000229        chrUn_gl000230 
#                  431                    67                   146 
#       chrUn_gl000231        chrUn_gl000232        chrUn_gl000233 
#                   92                   136                   154 
#       chrUn_gl000234        chrUn_gl000235        chrUn_gl000236 
#                  136                   115                   140 
#       chrUn_gl000237        chrUn_gl000238        chrUn_gl000239 
#                  153                   134                   113 
#       chrUn_gl000240        chrUn_gl000241        chrUn_gl000242 
#                  140                   141                   146 
#       chrUn_gl000243        chrUn_gl000244        chrUn_gl000245 
#                  145                   134                   123 
#       chrUn_gl000246        chrUn_gl000247        chrUn_gl000248 
#                  128                   122                   133 
#       chrUn_gl000249                  chrX                  chrY 
#                  129                517569                197912 



autosome = which(genometile[,1] %in% paste0("chr", 1:22))
length(autosome)
#[1] 9603454

autosome.tile = genometile[autosome,]
#path = "/mnt/gluster/data/external_private_supp/roger_atacseq/process_bams/"
#write.table(genometile[autosome,], file = paste0(path, "hg19.autosomes.tiles.w300.s300.bed"), quote=FALSE,row.names = FALSE, col.names = FALSE)
#write.table(autosome, file=paste0(path, "index.in.tiles.txt"), quote=FALSE, row.names=FALSE, col.names = FALSE)




max.posi = rep(NA, 22)
for(i in 1:22){
    wh = which(autosome.tile[,1] == paste0("chr",i))
    max.posi[i] = max(autosome.tile[wh,3])
}


max.posi
# [1] 249250621 243199373 198022430 191154276 180915260 171115067 159138663
# [8] 146364022 141213431 135534747 135006516 133851895 115169878 107349540
#[15] 102531392  90354753  81195210  78077248  59128983  63025520  48129895
#[22]  51304566

#path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC/info/"
#write.table(max.posi, file = paste0(path, "chr.len.txt"), quote=FALSE,row.names = FALSE, col.names = FALSE)





## Read logLR values from Tom
load("/mnt/lustre/home/shim/multiscale_analysis/test/etc/sim_res_ms_ss70.Robj")
loglr.ms.null.test[1:10]
## [1] 1.45225163 0.08760449 0.40487261 1.01319454 0.39671310 0.02202556
## [7] 1.16845057 0.92335796 0.27656655 1.84879880
loglr.ms.null[1:10]
## [1] 1.2429074 0.0783647 0.2915972 0.3191636 0.1379046 0.0000000 0.4757501
## [8] 0.1817970 0.2765666 1.5088881
loglr.ms.alt[1:10]
## [1]   7.9919357   0.0783647   2.1152773 129.8975699  23.3096496 455.0546313
## [7]  12.0434428  25.2813526  45.5932255  44.4825364


## Read logLR values from my analyses
path.null = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_manyQTLfinal_v1/null/multiscale/"
load(paste0(path.null, "sum/logLR.fullread.70ind.backup.Robj"))
logLR_list[1:10]
## [1] 0.332380678 0.078364696 0.290990494 0.293014182 0.003747977 0.000000000
## [7] 0.433140645 0.083485708 0.000000000 1.508888119

path.alt = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_manyQTLfinal_v2/alt/multiscale/"
load(paste0(path.alt, "sum/logLR.fullread.70ind.backup.Robj"))
logLR_list[1:10]
## [1] 11.5764757  0.0783647  4.8165593 34.3689072  0.4785586 58.5203679
## [7]  1.9593380  9.5663904  6.2289939 10.0315540


## Read logLR values from new analyses
path.null = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_manyQTLfinal_v1/null/multiscale/"
load(paste0(path.null, "sum/logLR.fullread.70ind.Robj"))
logLR_list[1:10]
## [1] 0.332380678 0.078364696 0.290990494 0.293014182 0.003747977 0.000000000
## [7] 0.433140645 0.083485708 0.000000000 1.508888119

path.alt = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_manyQTLfinal_v2/alt/multiscale/"
load(paste0(path.alt, "sum/logLR.fullread.70ind.Robj"))
logLR_list[1:10]
## [1] 11.5764757  0.0783647  4.8165593 34.3689072  0.4785586 58.5203679
## [7]  1.9593380  9.5663904  6.2289939 10.0315540

### Let's install a package using install.packages("~/multiseq/package/multiseq_0.1.tar.gz",repos=NULL,type="source")



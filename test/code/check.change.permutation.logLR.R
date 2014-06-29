
#### logLR computation 

numSites = 50


path.null = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/"

case.name = c("fullread.10ind","fullread.10ind.test", "fullread.10ind.test.pval")

logLR = vector("list", 3)


for(cc in 1:length(case.name)){

    logLR_list = rep(NA, numSites)
    done_res = rep(NA, numSites)

    IX = 1

    for(i in 1:numSites){
   
        num_path  = paste0(path.null, case.name[cc], ".output/res.", i, ".out")
        if(file.exists(num_path)== FALSE){		
            done_res[IX] = FALSE
        }else{
            if(file.info(num_path)$size == 0){
                done_res[IX] = FALSE
            }else{
                dat = scan(num_path, what="")
                if(!is.na(dat[1])){
                    logLR_list[IX] = as.numeric(dat[1])
                    done_res[IX] = TRUE
                }else{
                    done_res[IX] = FALSE
                }
            }
        }
        IX = IX + 1
    }

    logLR[[cc]] = logLR_list
    
    
}




sum(logLR[[1]] ==  logLR[[2]])
#[1] 50
sum(logLR[[2]] ==  logLR[[3]])
#[1] 50
sum(logLR[[3]] ==  logLR[[1]])
#[1] 50




### p-value computation 

numSites = 50


path.null = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/multiscale/"

case.name = c("fullread.10ind", "fullread.10ind.test.pval")

pval = vector("list", 2)


for(cc in 1:length(case.name)){

    pval_list = rep(NA, numSites)
    done_res = rep(NA, numSites)

    IX = 1

    for(i in 1:numSites){
   
        num_path  = paste0(path.null, case.name[cc], ".output/pval.", i, ".out")
        if(file.exists(num_path)== FALSE){		
            done_res[IX] = FALSE
        }else{
            if(file.info(num_path)$size == 0){
                done_res[IX] = FALSE
            }else{
                dat = scan(num_path, what="")
                if(!is.na(dat[4])){
                    pval_list[IX] = as.numeric(dat[4])
                    done_res[IX] = TRUE
                }else{
                    done_res[IX] = FALSE
                }
            }
        }
        IX = IX + 1
    }

    pval[[cc]] = pval_list
}

sum(pval[[1]] == pval[[2]])
#[1] 50



numSites = 578


path.alt = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/multiscale/"
path.null = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/null/multiscale/"

#case.name = c("fullread.70ind", "fullread.30ind", "fullread.10ind")
#case.name = c("fullread.70ind.over", "fullread.30ind.over", "fullread.10ind.over")
#case.name = c("fullread.70ind.over.2", "fullread.30ind.over.2", "fullread.10ind.over.2")
#case.name = c("halfread.70ind.over", "halfread.30ind.over")
#case.name = c("2fullread.70ind.over", "2fullread.30ind.over", "2fullread.10ind.over", "4fullread.10ind.over")
case.name = c("4fullread.70ind.over", "4fullread.30ind.over")



done_list_alt = vector("list", length(case.name))
done_list_null = vector("list", length(case.name))


for(cc in 1:length(case.name)){

    pval_list = rep(NA, numSites)
    done_res = rep(NA, numSites)

    IX = 1

    for(i in 1:numSites){
   
        num_path  = paste0(path.alt, case.name[cc], ".output/pval.", i, ".out")
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
    done_list_alt[[cc]] = done_res
    save("pval_list", "done_res", file =paste0(path.alt, "sum/pval.", case.name[cc], ".Robj"))


}




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
    done_list_null[[cc]] = done_res
    save("pval_list", "done_res", file =paste0(path.null, "sum/pval.", case.name[cc], ".Robj"))

    
}





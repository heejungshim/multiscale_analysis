

numSites = 500


path.alt = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/alt/wave/"
path.null = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/null/wave/"


case.name = c("halfread.70ind.over", "halfread.30ind.over", "2fullread.10ind.over")


#case.name = c("fullread.70ind.over", "fullread.30ind.over", "fullread.10ind.over")
#case.name = c("fullread.70ind.over.2", "fullread.30ind.over.2", "fullread.10ind.over.2")
#case.name = c("fullread.70ind", "fullread.30ind", "fullread.10ind")



for(cc in 1:length(case.name)){

    #cc = 1
    logLR_list = rep(NA, numSites)
    done_res = rep(NA, numSites)

    IX = 1

    num_dir_path = paste0(path.alt, "output/", case.name[cc], ".")

    for(i in 1:numSites){
        # i = 1
        num_path  = paste0(num_dir_path, i, ".fph.BFs.txt")
        if(file.exists(num_path)== FALSE){		
            done_res[IX] = FALSE
        }else{
            if(file.info(num_path)$size == 0){
		done_res[IX] = FALSE
            }else{			
		dat = read.table(num_path, as.is = TRUE)
                if(!is.na(dat[1,2])){
                    logLR_list[IX] = as.numeric(dat[1,2])
                    done_res[IX] = TRUE
                }else{
                    done_res[IX] = FALSE
                }
            }
        }
        IX = IX + 1
    }

    save("logLR_list", "done_res", file =paste0(path.alt, "sum/logLR.", case.name[cc], ".Robj"))

}



for(cc in 1:length(case.name)){

    logLR_list = rep(NA, numSites)
    done_res = rep(NA, numSites)

    IX = 1

    num_dir_path = paste0(path.null, "output/", case.name[cc], ".")

    for(i in 1:numSites){
   
        num_path  = paste0(num_dir_path, i, ".fph.BFs.txt")
        if(file.exists(num_path)== FALSE){		
            done_res[IX] = FALSE
        }else{
            if(file.info(num_path)$size == 0){
		done_res[IX] = FALSE
            }else{			
		dat = read.table(num_path, as.is = TRUE)
		if(!is.na(dat[1,2])){
                    logLR_list[IX] = as.numeric(dat[1,2])
                    done_res[IX] = TRUE
                }else{
                    done_res[IX] = FALSE
                }
            }
        }
        IX = IX + 1
    }

    save("logLR_list", "done_res", file =paste0(path.null, "sum/logLR.", case.name[cc], ".Robj"))

}



numSites = 578


path.alt = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/wave/"
path.null = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/null/wave/"


#case.name = c("halfread.70ind.over", "halfread.30ind.over")


#case.name = c("fullread.70ind.over", "fullread.30ind.over", "fullread.10ind.over", "fullread.4ind.over")
#case.name = c("fullread.70ind.over.2", "fullread.30ind.over.2", "fullread.10ind.over.2", "fullread.4ind.over.2")
#case.name = c("fullread.70ind", "fullread.30ind", "fullread.10ind", "fullread.4ind")

#case.name = c("2fullread.4ind.over", "4fullread.4ind.over")
#case.name = c("fullread.6ind.over", "2fullread.6ind.over", "4fullread.6ind.over")

#case.name = c("2fullread.70ind.over", "2fullread.30ind.over", "2fullread.10ind.over", "4fullread.10ind.over")

case.name = c("4fullread.70ind.over", "4fullread.30ind.over")



done_list_alt = vector("list", length(case.name))
done_list_null = vector("list", length(case.name))


for(cc in 1:length(case.name)){

    #cc = 1
    logLR_list = rep(NA, numSites)
    done_res = rep(NA, numSites)

    IX = 1

    num_dir_path = paste0(path.alt, "output/", case.name[cc], ".")

    for(i in 1:numSites){
        # i = 1
        num_path  = paste0(num_dir_path, i, ".fph.logLR.txt") 
        #num_path  = paste0(num_dir_path, i, ".fph.BFs.txt")
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
    done_list_alt[[cc]] = done_res
    save("logLR_list", "done_res", file =paste0(path.alt, "sum/logLR.", case.name[cc], ".Robj"))

}



for(cc in 1:length(case.name)){

    logLR_list = rep(NA, numSites)
    done_res = rep(NA, numSites)

    IX = 1

    num_dir_path = paste0(path.null, "output/", case.name[cc], ".")

    for(i in 1:numSites){
   
        #num_path  = paste0(num_dir_path, i, ".fph.BFs.txt")
        num_path  = paste0(num_dir_path, i, ".fph.logLR.txt")

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
    done_list_null[[cc]] = done_res
    save("logLR_list", "done_res", file =paste0(path.null, "sum/logLR.", case.name[cc], ".Robj"))

}

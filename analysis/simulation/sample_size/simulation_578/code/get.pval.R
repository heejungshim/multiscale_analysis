

numSites = 578
#path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/alt/multiscale/"
path = "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_578/null/multiscale/"



pval_list = rep(NA, numSites)
done_res = rep(NA, numSites)

IX = 1

for(i in 1:numSites){
   
    num_path  = paste0(path, "output/pval.", i, ".out")
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


sum(done_res)
max(as.numeric(pval_list))
min(as.numeric(pval_list))
which(done_res == 0)
#[1]  12  86 334 431 507 525 562
save("pval_list", "done_res", file =paste0(path, "sum/pval.Robj"))







In simulation with effect size from multiseq (see compareWaveQTLandmultiseq_manydsQTL_final.org): case with `halfread.10ind' and `seed=225'
See the scipt `simulateDataForAll.run.multiscale_debug.R' in ~/multiscale_analysis/test/code/

phenoD : phenotype data for 10 individuls on 1024bp
genoD : group indicator
res : result from multiseq (multiseq version: 0.1.0 and ashr version: 0.9.3)
command : res=multiseq(phenoD,genoD,ashparam=list(prior="uniform"),verbose=TRUE)

Issue is logLR from 10th scale is NA. 


> load("data.Robj")
> ls()
 [1] "genoD"  "phenoD" "res" 
> dim(phenoD)
 [1]   10 1024
> genoD
 [1] 0 0 0 0 0 1 1 1 1 1
> str(res)
List of 7
 $ baseline.mean     : num [1:1024] NA NA NA NA NA NA NA NA NA NA ...
 $ baseline.var      : num [1:1024] NA NA NA NA NA NA NA NA NA NA ...
 $ effect.mean       : num [1:1024] NA NA NA NA NA NA NA NA NA NA ...
 $ effect.var        : num [1:1024] NA NA NA NA NA NA NA NA NA NA ...
 $ logLR             :List of 3
  ..$ value   : num NA
  ..$ scales  : num [1:11] 0 0 0 0 0 0 0 0 0 NA ...
  ..$ isfinite: logi FALSE
 $ fitted.g          :List of 11
  ..$ : NULL
  ..$ : NULL
  ..$ : NULL
  ..$ : NULL
  ..$ : NULL
  ..$ : NULL
  ..$ : NULL
  ..$ : NULL
  ..$ : NULL
  ..$ : NULL
  ..$ :List of 3
  .. ..$ pi  : num [1:5] 1 0 0 0 0
  .. ..$ mean: num [1:5] 0 0 0 0 0
  .. ..$ sd  : num [1:5] 0 695 1390 2780 5560
  .. ..- attr(*, "row.names")= int [1:5] 1 2 3 4 5
  .. ..- attr(*, "class")= chr "normalmix"
 $ fitted.g.intercept: list()
 - attr(*, "class")= chr "multiseq"



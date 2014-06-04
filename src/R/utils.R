## This file contains useful functions I copied from somewhere (see each function for sources)




## 'estBetaParams' takes mean and variance for beta distribution and returns alpha and beta parameters.  
## I copied the original version from http://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
## Then I modified to handle cases when somme elements of mu are 0 or 1 (then it's always 0 or 1 as I don't allow variation.)
estBetaParams <- function(mu.orig, var) {

    del.ix = ((mu.orig == 0) | (mu.orig == 1))
    if(sum(del.ix) > 0){
        mu = mu.orig[!del.ix]
    }else{
        mu = mu.orig
    }
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)

    if(sum(del.ix) > 0){
        alpha.orig = beta.orig = rep(NA, length(mu.orig))
        alpha.orig[!del.ix] = alpha
        beta.orig[!del.ix] = beta
    }else{
        alpha.orig = alpha
        beta.orig = beta
    }
    return(params = list(alpha = alpha.orig, beta = beta.orig))
}
















## This file contains useful functions I copied from somewhere (see each function for sources)




## 'estBetaParams' takes mean and variance for beta distribution and returns alpha and beta parameters.  
## I copied from http://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
estBetaParams <- function(mu, var) {
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))
}






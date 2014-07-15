## This file contains useful functions for analysis.
## Sometime I copied original functions and modified them. Then, I specified original sources. 
##
## Copyright (C) 2014 Heejung Shim
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.




##' 'estBetaParams' takes mean and variance for beta distribution and returns
##' alpha and beta parameters in beta distribution.
##'
##'
##' I copied the original version from http://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
##' Then I modified to handle cases when somme elements of mu are 0 or 1 (then it returns NA if mu <= 0 or >=1) or when computed alpha and beta <= 0 (then it returns NA)
##'
##' @param mu.orig a vector of mean for beta distribution
##' @param var a vector (or scalar) of variance for beta distribution
##' @return a list of alpha and beta 
estBetaParams <- function(mu.orig, var) {

    del.ix = ((mu.orig <= 0) | (mu.orig >= 1))
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

    ## handle non-positive alpha or beta 
    invalid.para = which((alpha.orig <= 0) | (beta.orig <= 0))
    if(length(invalid.para) > 0){
        alpha.orig[invalid.para] = rep(NA, length(invalid.para))
        beta.orig[invalid.para] = rep(NA, length(invalid.para))
    }
    
    return(params = list(alpha = alpha.orig, beta = beta.orig))
}





##' 'sample.from.Binomial.with.Overdispersion' simulates binomial samples with/without
##' over dispersion. 
##'
##' For a given overdispersion parameter, computed parameters for beta distribution can be invalid (e.g., mu.sig are 0 or 1). Then we sample read without overdispersion for those positions.
##' 
##' @param num.sam number of samples to be sampled
##' @param total.count a vector of non-negative counts;
##' @param mu.sig a vector of probabilities (we allow 0 or 1 as probablity)
##' @param over.dispersion if over.dispersion == NULL, simulate data from binomial. If over.dispersion is provided, simulate data binomial with over.dispersion.
##' @return a matrix of num.sam by L (length of total.count) containing simulated data. 
sample.from.Binomial.with.Overdispersion <- function(num.sam, total.count, mu.sig, over.dispersion=NULL){

    invalid.entry = ((mu.sig < 0) | (mu.sig > 1))
    if(sum(invalid.entry) > 0){ stop("ERROR, mu.sig have some values outside of valid range [0, 1]")}
                                
     if(is.null(over.dispersion)){
        return(matrix(data=rbinom(length(mu.sig)*num.sam, total.count, mu.sig), nr = num.sam, byrow = TRUE))
    }else{

        
        final.dat = matrix(data=NA, nr = num.sam, nc = length(mu.sig))
        
        # get alpha and beta
        resBeta = estBetaParams(mu.sig, over.dispersion)
        alpha = resBeta$alpha
        beta = resBeta$beta

        # for valid alpha and beta, sample data 
        del.ix = is.na(alpha)
        p.sig = rbeta(sum(!del.ix)*num.sam, alpha[!del.ix], beta[!del.ix])
        dat.V = rbinom(sum(!del.ix)*num.sam, total.count[!del.ix], p.sig) 
        final.dat[,!del.ix] = matrix(data=dat.V, nr = num.sam, byrow = TRUE)

        # for invalid alpha and beta, sample without over dispersion
        if(sum(del.ix) > 0){
            dat.IV = matrix(data=rbinom(sum(del.ix)*num.sam, total.count[del.ix], mu.sig[del.ix]), nr = num.sam, byrow = TRUE)
            final.dat[,del.ix] = matrix(data=dat.IV, nr = num.sam, byrow = TRUE)
        }

        return(final.dat)

    }

}
 













# function to compute alpha and beta from mean and varince of beta distribution
# copyed from http://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
estBetaParams <- function(mu, var) {
     alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
     beta <- alpha * (1 / mu - 1)
     return(params = list(alpha = alpha, beta = beta))
}


mu = 1/70
var = 1/70/70/2

alpha = estBetaParams(mu, var)$alpha
beta = estBetaParams(mu, var)$beta

sam1 = rbeta(10000, alpha, beta)


mu = 1/70
var = 1/70/70/4

alpha = estBetaParams(mu, var)$alpha
beta = estBetaParams(mu, var)$beta

sam2 = rbeta(10000, alpha, beta)



mu = 1/70
var = 1/70/70/7

alpha = estBetaParams(mu, var)$alpha
beta = estBetaParams(mu, var)$beta

sam3 = rbeta(10000, alpha, beta)


## make a histrogam

pdf("beta.samples_v2.pdf")
par(mfrow=c(3,1))
hist(sam1, breaks =100, main = "mean: 1/70, var: 1/70/70/2")
abline(v = 1/70, col="red")

hist(sam2, breaks =100, main = "mean: 1/70, var: 1/70/70/4")
abline(v = 1/70, col="red")

hist(sam3, breaks =100, main = "mean: 1/70, var: 1/70/70/7")
abline(v = 1/70, col="red")

dev.off()



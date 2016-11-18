#########################################################
## Bayesian Statistics: MCMC examples
## Instructor: Jason Guo
## Quant III Lab 10
#########################################################

##################
## Gibbs Sampler #
##################

y <- rnorm(20, 5, 20)

gibbs <- function(n.sim, mu0, tau2, a, b){
  n <- length(y)
  theta     <- list(mu=0, sigma2=1)      # Initial values
  theta.sim <- matrix(NA, n.sim, 2)
  for (i in 1:n.sim){
    theta$mu     <- rnorm(1, (mu0/tau2 + n*mean(y)/theta$sigma2)/(1/tau2 + n/theta$sigma2), sqrt(1/(1/tau2 + n/theta$sigma2)))
    theta$sigma2 <- 1/rgamma(1, n/2 + a , .5*sum((y-theta$mu)^2) + b) 
    theta.sim[i, ]  <- unlist(theta)
  }
  theta.sim
}

# posterior distribution
plot(density(gibbs(1000, 100, 100, 5, 5)[,1]))
# posterior mean
mean(gibbs(1000, 100, 100, 5, 5)[,1])
plot(density(gibbs(1000, 30, 100, 5, 5)[,2]))
#########################################################
## Bayesian Statistics: MCMC examples
## Instructor: Jason Guo
## Quant III Lab 10
#########################################################

################
# markov chain #
################

#Every person in a society can be either sick (S), healthy but not immune (H), or immune
#and healthy (I). People move from one state to another according to the following Markov
#transition probabilities:

transition <- matrix(c(0.9,0,0.1,0.6,0.4,0,0,0,1),nrow=3,byrow=TRUE)

# Plot the distribution of the states as a function of n, starting with a state where 50% of
# people are healthy and 50% are sick.

p <- as.vector(c(0, 0.5,0.5))

x <- rep(NA, 100)
y <- rep(NA, 100)
z <- rep(NA, 100)

x[1] <- 0
y[1] <- 0.5
z[1] <- 0.5

for (i in 1:99){ 
  d <- t(p) %*% transition # distribution at t = i + 1, using transition matrix at t = i
  x[i+1] <- d[1]
  y[i+1] <- d[2]
  z[i+1] <- d[3]
  transition <- transition %*% transition # transition matrix at t = i + 1
}

#distribution <- cbind(x,y,z)
#plot.ts(distribution)

plot(x, type="l", col=1, ylim=c(0,1), lwd=2,
     ylab="Probability", main="Distribution Across Time",xlab="Time")
lines(y, col=2,lwd=2 ,lty=2)
lines(z, col=4, lwd=2)
legend("right", legend=c("Sick","Healthy","Immune"), col=c(1,2,4), lty=c(1,2,1), lwd=2)

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
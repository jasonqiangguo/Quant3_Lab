#########################################################
## Bayesian Statistics: MCMC examples
## Instructor: Jason Guo
## Quant III Lab 10
#########################################################

library(pscl)

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
# classical MCMC #
##################

##########################################################################################
# the prior for \mu is a normal distribution with mean mu0 and variance tau2, the prior  #
# for \sigma2 is an inverse-gamma distribution with parameter (a, b), and the likelihood #
# is a normal distribution with mean \mu and variance \sigma2                            #
##########################################################################################

##################
## Gibbs Sampler #
##################


y <- rnorm(100, 5, 5)

gibbs <- function(mu0, tau2, a, b, n.sim){
  n <- length(y)
  theta     <- list(mu=0, sigma2=1)      # Initial values
  theta.sim <- matrix(NA, n.sim, 2)
  for (i in 1:n.sim){
    theta$mu     <- rnorm(1, (mu0*theta$sigma2 + n*mean(y)*tau2)/(n*tau2 + theta$sigma2), sqrt(1/(1/tau2 + n/theta$sigma2)))
    theta$sigma2 <- rigamma(1, n/2 + a , .5*sum((y-theta$mu)^2) + b) 
    theta.sim[i, ]  <- unlist(theta)
  }
  theta.sim
}


########################
## metropolis hastings #
########################

# likelihood of joint posterior p(\mu, \sigma2|y) \propto p(y|\mu, \sigma2)p(\mu)p(\sigma2)

# one caveat is that \sigma2 is positive, so in order to avoid negative draws for \sigma2, we instead draw t
# such that exp(t) = \sigma2. In addition, to avoid numeric overflow, we take the log of the likelihood for
# the posterior

# loglikelihood of the joint posterior 
lpost <- function(theta, mu0, tau2, a, b){
  log.post <- -(length(y/2) + a + 1)*theta[2] - (theta[1] - mu0)^2/(2*tau2) - (2*b + sum((y- theta[1])^2))/(2*exp(theta[2]))
}

# random walk m-h algorithm
mh <- function(theta, mu0, tau2, a, b ,s) {
  cand <-  theta + rnorm(length(theta), 0, s)
  accept <- lpost(cand, mu0, tau2, a, b) - lpost(theta, mu0, tau2, a, b)
  accept <- exp(min(accept, 0))
  if(runif(1) <= accept) theta <- cand
  theta
}

# burnin is the period of instability before getting into convergence 
# thin is the parameter specifying the distance between each draw to be stored in the matrix
metro <- function(mu0, tau2, a, b, s, n.sim, burnin, thin){
  theta     <- rnorm(2) 
  theta.sim <- matrix(NA, n.sim, 2)  # matrix to store values
  for (i in 1:n.sim){
    theta          <- mh(theta, mu0, tau2, a, b , s)
    theta.sim[i, ] <- c(theta[1], exp(theta[2]))
  }
  thin.draw <- seq((burnin + 1), n.sim, by = thin)
  theta.sim[thin.draw,]
}

ytry <- rnorm(100, 0, 1)
xtry <- rnorm(100, 2, 3)
m <- as.matrix(cbind(ytry, xtry))
h <- seq(4, 80, by = 3)
m[h,]

m.data <- data.frame(metro = metro(-5, 1, 2, 2, 2, n.sim = 1000000, burnin = 20000, thin = 4))
g.data <- data.frame(gibbs = gibbs(-5, 1, 2, 2, 10000))

plot(density(m.data$metro.1))
plot(density(m.data$metro.2))

plot(density(pdata$gibbs.1))
plot(density(pdata$gibbs.2))


#################################################################################################


#########
# Rstan #
#########

library(rstan)

hmc_code <- '
data {
int<lower=0> N;
vector[N] y ;
}

parameters {
real mu;
real<lower=0> sigma2; 
}

model {
mu ~ normal(-5, 1);
sigma2 ~ inv_gamma(2,2);
for (n in 1:N)
y[n] ~ normal(mu, sqrt(sigma2));
}
'

data <- list(N=length(y), y=y)
fit <- stan(model_code=hmc_code, data=data, iter=5000, warmup = 200, thin = 3, chains=4)
monitor(fit)
traceplot(fit, pars='mu')

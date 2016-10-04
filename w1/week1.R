################################################
## Estimating Possion Distribution Parameters 
## Instructor: Jason Guo
## Quant III Lab 1
## Sep 15th
################################################

library(foreign)
library(dplyr)

rm(list=ls())
setwd("/Users/qiangguo/Google Drive/Quant3_TA/TA/w1")

# First, let us do a small exercise to recover parameters of poisson distribution with MLE
# Data Generating Process Y ~ Poisson (lambda)
Y <- rpois(500, 3)
summary(Y)
plot(table(Y))

# Write the joint likelihood function (a function of the parameter lambda) to be maximized
poisson.lk <- function(lambda){
  l <- sum(lambda - Y * log(lambda) + log(factorial(Y)))
  return(l)
}

# or

poisson.lk <- function(lambda){
  l <- sum(lambda - Y * log(lambda))
  return(l)
}

# These two likelihood functions give the same results as the sum of log(factorial(Y)) is fixed and will not affect the maximization(minimization)

# Use optim to maximize
ml.estimates <- optim(par=1, poisson.lk, method="BFGS", hessian = TRUE)
# check parameters
ml.estimates$par
# check Hessian
ml.estimates$hessian
# obtain covariance-variance matrix
cov <- solve(ml.estimates$hessian)
# obtain standard errors
lambda.se <- sqrt(diag(cov))


# Another way to construct likelihood funciton is that we sum the logged density/probability given we know exactly the distribution
# Remember the object we are to maximize is the joint likelihood of the parameter given all the data points, which is the product of
# all the probabilities across the span of data. 

poisson.lk <- function(lambda){
  l <- -sum(dpois(Y, lambda, log = TRUE))
  return(l)
}

# to see whether dpois gives the likelihood
plot(Y, dpois(Y, 3))


# assessing convergence of ML estimator using simulation
n <- seq(5, 5000, 25)

# rewrite poisson joint likelihood function
pois.lk <- function(lambda, y){
  l <- -sum(dpois(y, lambda, log = TRUE))
  return(l)
}

convergence <- function(n){
  y <- rpois(n, 5)
  ml.estimates <- optim(par=1, pois.lk, y = y, method="BFGS")
  return(ml.estimates$par)
}

results <- sapply(n, convergence)

plot(n, results, type="l", ylim=c(4, 6), 
     xlab="Sample size", ylab="Estimated lambda")
abline(h=5, col="grey80")

plot(n, results-5, type="l", ylim=c(-1, 1), 
     xlab="Sample size", ylab="lambda_hat - lambda")
abline(h=0, col="grey80")

# visualizing log.lk distribution
y <- rpois(100, 5)
lambdas <- seq(1, 10, 0.10)
lks <- -sapply(lambdas, pois.lk, y)

plot(lambdas, lks, type="l", xlab="lambda", ylab="logL")




# Now, let us use real data to estimate parameters with MLE
# Fearon and Latin (2003)
# They want to explore the driving factors of occurences of conlicts in a country. Here we assume that ethnic fractionalization, population and their interaction may have
# an effect on the number of conlicts/wars a country experiences.

d1 <- read.dta("conflict.dta")
# There are five variables in the dataset, ccode (country code), war (number of wars a country experienced), ethfrac (measure of average ethnic diversity),
# lpop (average logged population), interaction of averaged ethfrac and lpop, and constant term


# Plotting the distribution of wars, you will find excessive zeros, which is a problem as poisson distribution will not produce so disproportionately many zeros.
# But here we temporarily ignore it and we will come back to this issue in future sessions.
plot(table(d1$war))



set.seed(333)
# set initials for ML, here we have 4 parameters to be estimated, therefore we need to have 4 initial values. 
initial <- rnorm(4,mean=0,sd=1)

log.likelihood <- function(betas){
  lambda <- exp(d1$constant*betas[1] + d1$ethfrac*betas[2] + d1$lpop*betas[3] + d1$interaction*betas[4])
  l <- -sum(dpois(d1$war, lambda, log = TRUE))
}

mle <- optim(initial, fn=log.likelihood, method="BFGS",hessian=TRUE)
mleest <- mle$par
cov <- solve(mle$hessian)
betase <- sqrt(diag(cov))

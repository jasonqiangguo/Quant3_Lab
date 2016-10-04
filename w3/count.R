rm(list=ls())
library(foreign)

# Overdispersion and zero inflation - POISSION and BINOMIAL
data <- read.dta("soccer_data.dta")
# Miguel, Edward, Sebastián M. Saiegh, and Shanker Satyanath. 
# "Civil war exposure and violence." Economics & Politics 23.1 (2011): 59-73.

## estimating poisson model
#---------------------------------------------------------------------
summary(
  poisson <- glm(yellow_card ~ civwar + age, family = "poisson", data = data)
)

## negative binomial
library(MASS)
summary(
  negbin <- glm.nb(yellow_card ~ civwar + age, data = data)
)
## estimating marginal effects of age using simulation, assuming civwar=0
#---------------------------------------------------------------------
par(mfrow=c(1,2))
## POISSON
betas <- mvrnorm(n=1000, mu=poisson$coefficients, 
                 Sigma=summary(poisson)$cov.unscaled)
age <- 18:40
marg <- matrix(NA, nrow=1000, ncol=length(age))
for (i in 1:1000){
  for (a in 1:length(age)){
    marg[i,a] <- betas[i,3] * exp(c(1, 0, a) %*% betas[i,])
  }
}

marg <- apply(marg, 2, function(x)
  c(mean(x), quantile(x, c(.025, .975))))

plot(age, marg[1,], type="l", 
     xlab="Age", ylab="Marginal effect of age", ylim=c(0.02, 0.07))
lines(age, marg[2,], col="grey")
lines(age, marg[3,], col="grey")

#---------------------------------------------------------------------
## NEGATIVE BINOMIAL
betas <- mvrnorm(n=1000, mu=negbin$coefficients, 
                 Sigma=summary(negbin)$cov.unscaled)
age <- 18:40
marg <- matrix(NA, nrow=1000, ncol=length(age))
for (i in 1:1000){
  for (a in 1:length(age)){
    marg[i,a] <- betas[i,3] * exp(c(1, 0, a) %*% betas[i,])
  }
}

marg <- apply(marg, 2, function(x)
  c(mean(x), quantile(x, c(.025, .975))))

plot(age, marg[1,], type="l", 
     xlab="Age", ylab="Marginal effect of age", ylim=c(0.02, 0.07))
lines(age, marg[2,], col="grey")
lines(age, marg[3,], col="grey")
#---------------------------------------------------------------------
# Infalted Zeros?
par(mfrow=c(1,1))
hist(data$yellow_card)

# How would you test for it?











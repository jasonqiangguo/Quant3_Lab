#############################################################
## Bayesian Statistics: Multilevel model/ Hierarchical Model
## Instructor: Jason Guo
## Quant III Lab 11
#############################################################

rm(list=ls())


library(pscl)
library(rstan)
library(coda)
library(ggplot2)
library(ggmcmc)

# Stein's Paradox

# this example is based on Jackman 310-316 and Pablo's code

data(EfronMorris)
?EfronMorris
d <- EfronMorris


# MLE estimate is just batting average
rmse <- sqrt(mean((d$y-d$p)^2))
rmse


# Hierarchical model ('partial pooling')
#---------------------------------------------

# each observed batting average y_i is drawn from a normal density with
# unknown mean theta_i; 


# y_i~N(\theta_i,sigma^2)
# \theta_i~N(\mu,\tau^2)
# \tau^2~IG(a,b)
# \mu, a,b are per-specified hyperparameter

# and fixed variance sigma:
sigma <- sqrt(mean(d$y) * (1-(mean(d$y)))/45)
sigma


# prior on \mu: we know it must be somewhere between 0.15 and 
# .30 (read the example); so we choose values that give 
# a 95% credible interval: normal(0.225, 0.0375)

# what prior belief do we have about tau:
# the standard deviation of the thetas?
# the sd of the observed y is:
sd(d$y)

# it must be somewhere around 0.07; probably higher; let's say 0.10
# so MEAN of values we draw from gamma distribution must be around 0.10
# gamma distribution has mean shape/rate [OR shape*scale]

# so let's follow Jackman and choose shape=7; then rate=7/0.10 ~= 70
x <- seq(0.001, 1, 0.001)
plot(x, dgamma(x, shape=7, rate=70), type="l")


hierarchical_code <- '
data {
int<lower=0> N; # observations
real y[N]; # outcome variable
real<lower=0> sigma; # fixed variance of observed batting averages
}
parameters {
real theta[N]; ## unobserved true batting averages
real mu; ## hyperparameter: mean of batting averages
real<lower=0> tau; ## hyperparameter: variance of batting averages
}
model {
mu ~ normal(0.225, 0.0375); ## prior about hyperparameter (overall batting average)
tau ~ gamma(7, 70); ## priors about variance of overall batting average
for (n in 1:N){;
y[n] ~ normal(theta[n], sigma);
theta[n] ~ normal(mu, tau);
}
}
'

data <- list(N=length(d$y), y=d$y, sigma=sigma)

n.chain <- 4
myinits <- list()
for (i in 1:n.chain) {
  myinits[[i]] <- list(theta = rep(0, 18), mu = 0, tau = 0.1)
}

fit <- stan(model_code=hierarchical_code, data=data, iter=10000, chains=4, init = myinits)


# summary results
monitor(fit, digits_summary=3)

# Assessing convergence
mcmc <- As.mcmc.list(fit, pars = "theta")
S <- ggs(mcmc, inc_warmup = TRUE)
theta_1_to_9 <- dplyr::filter(S, Parameter == "theta[1]"|Parameter == "theta[2]"|Parameter == "theta[3]"|
                                Parameter == "theta[4]"|Parameter == "theta[5]"|Parameter == "theta[6]"|
                                Parameter == "theta[7]"|Parameter == "theta[8]"|Parameter == "theta[9]")

ggs_running(theta_1_to_9)

ggs_traceplot(theta_1_to_9, original_burnin = TRUE)

ggs_autocorrelation(theta_1_to_9)

ggs_compare_partial(theta_1_to_9)

ggs_geweke(S) #between -2 and 2 is good# 


# Heidelberger-Welch test of non-stationarity
heidel.diag(mcmc)


# computing rmse
estimates <- summary(fit)
predicted <- estimates$summary[1:18, 'mean']
rmse <- sqrt(mean((predicted-d$p)^2))
rmse

# MLE estimate is just batting average
rmse <- sqrt(mean((d$y-d$p)^2))
rmse


# visualizing results (replication of Figure 7.4 in page 316 Jackman)
actual <- data.frame(batter = d$name, average=d$p, estimate="Actual")
bayes <- data.frame(batter = d$name, average=predicted, estimate="Bayes")
mle <- data.frame(batter = d$name, average=d$y, estimate="MLE")

df <- rbind(actual, bayes, mle)

p <- ggplot(df, aes(x=estimate, y=average, group=batter))
p + geom_point() + geom_line() + theme_bw() + coord_flip() +
  scale_x_discrete(expand=c(0.05,0.05)) +
  scale_y_continuous("Batting average") + 
  theme(axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size = 0.50)) +
  geom_vline(xintercept=1, alpha=1/3, size=0.2) +
  geom_vline(xintercept=2, alpha=1/3, size=0.2) +
  geom_vline(xintercept=3, alpha=1/3, size=0.2)







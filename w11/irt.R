#########################################################
## Bayesian Statistics: IRT model
## Instructor: Jason Guo
## Quant III Lab 11
#########################################################

rm(list=ls())

library(pscl)
library(rstan)
library(ggplot2)
library(ggstance)
library(ggmcmc)
library(coda)
set.seed(2016)


setwd("~/GitHub/Quant3_Lab/w11")
# loading roll call data
# Source: http://jackman.stanford.edu/blog/
load("rc.rda")
rc <- dropRollCall(rc, dropList=list(codes = "notInLegis", lop = 0))


stan.code <- '
data {
  int<lower=1> J; // number of legislators
  int<lower=1> K; // number of bills
  int<lower=1> N; // number of observations
  int<lower=1,upper=J> j[N]; // legislator for observation n
  int<lower=1,upper=K> k[N]; // bill for observation n
  int<lower=0,upper=1> y[N]; // vote of observation n
}
parameters {
  real alpha[K];             
  real beta[K];   
  real theta[J];
}
model {
  alpha ~ normal(0, 25);
  beta ~ normal(0, 25);
  theta ~ normal(0, 1);
  for (n in 1:N)
    y[n] ~  bernoulli_logit(theta[j[n]] * beta[k[n]] - alpha[k[n]]);
}
'

J <- dim(rc$votes)[1]
K <- dim(rc$votes)[2]
N <- length(rc$votes)
j <- rep(1:J, times=K)
k <- rep(1:K, each=J) 
y <- c(rc$votes)

# deleting missing values
miss <- which(is.na(y))
N <- N - length(miss)
j <- j[-miss]
k <- k[-miss]
y <- y[-miss]

## data and initial values
stan.data <- list(J=J, K=K, N=N, j=j, k=k, y=y)

n.chain <- 3
myinits <- list()
for (i in 1:n.chain) {
  myinits[[i]] <- list(alpha = rep(0.5, K), beta = rep(0.5, K), theta = rep(0.5, J))
}

stan.fit <- stan(model_code=stan.code, data=stan.data, iter=1000, warmup=200,
                 chains= 3, thin=2, init = myinits)


##########################################################################

load("stan_fit_irt_1000.Rdata")


# convergence check, use "coda" and "ggmcmc"

# convert a stan object to a coda object
s <- As.mcmc.list(stan.fit, pars = c("theta"))

# use ggs to convert coda object to ggmcmc object for graphics
S <- ggs(s, inc_warmup = TRUE)
theta14 <- dplyr::filter(S, Parameter == "theta[14]")
ggs_running(theta14)
ggs_traceplot(theta14, original_burnin = TRUE)
ggs_compare_partial(theta14)
ggs_autocorrelation(theta14)
ggs_geweke(S)

#Heidelberger-Welch test of non-stationarity
heidel.diag(S)

# plot the ideological positions of senators
HPD <- data.frame(HPDinterval(mcmc(rbind(s[[1]], s[[2]], s[[3]]))))
plot.data <- as.data.frame(monitor(stan.fit))
plot.data.theta <- plot.data[grep("theta", rownames(plot.data)),]
plot.data.theta$name <- rownames(rc$votes)
plot.data.theta$partyid <- rc$legis.data$party
plot.data.theta$lower <- HPD$lower
plot.data.theta$upper <- HPD$upper

p <- ggplot(plot.data.theta, aes(y= reorder(name, -mean), x = mean)) + geom_point() +
  geom_linerangeh(aes(xmin= lower, xmax= upper, color=partyid, linetype = partyid)) + 
  theme_bw() + theme(axis.title.y=element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), legend.position = c(0.15,0.2),
                     plot.title = element_text(hjust = 0.4, size = rel(2)), legend.text = element_text(size = 15),
                     legend.title = element_text(size = 20)) +
  scale_color_manual(name = "Party ID", values = c("D" = "blue", "I" = "black", "R" = "red"), labels = c("Democrat", "Independent", "Republican")) + 
  scale_linetype_manual(name = "Party ID", values = c("D" = 1, "I" = 2, "R" = 4), labels = c("Democrat", "Independent", "Republican")) + 
  labs(x = "Ideological Position with 95% HPD Interval", title = "113th Congress: Senate")  
p
ggsave("113senate.pdf", height = 12, width = 8)


##############################################################################################################
# if you find democrats on the right of spectrum, you can reverse data. This issue can rise due to reflection
# problem, as theta_i * b_j - d_j = -theta_i * -b_j - d_j clearly leads to multimode for theta and b_j
##############################################################################################################

reverse.plot.data.theta <- plot.data.theta
reverse.plot.data.theta$mean <- -plot.data.theta$mean
reverse.plot.data.theta$lower <- -plot.data.theta$lower
reverse.plot.data.theta$upper <- -plot.data.theta$upper

p <- ggplot(reverse.plot.data.theta, aes(y= reorder(name, -mean), x = mean)) + geom_point() +
  geom_linerangeh(aes(xmin= lower, xmax= upper, color=partyid, linetype = partyid)) + 
  theme_bw() + theme(axis.title.y=element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), legend.position = c(0.15,0.2),
                     plot.title = element_text(hjust = 0.4, size = rel(2)), legend.text = element_text(size = 15),
                     legend.title = element_text(size = 20)) +
  scale_color_manual(name = "Party ID", values = c("D" = "blue", "I" = "black", "R" = "red"), labels = c("Democrat", "Independent", "Republican")) + 
  scale_linetype_manual(name = "Party ID", values = c("D" = 1, "I" = 2, "R" = 4), labels = c("Democrat", "Independent", "Republican")) + 
  labs(x = "Ideological Position with 95% HPD Interval", title = "113th Congress: Senate")  
p
##########################################################################################################
ggsave("113senate.pdf", width = 8, height = 12)


################################################################
## IRT WITH COVARIATES
################################################################

# Different ways of doing this...
# With Stan, it would be something like this:


stan.code <- '
data {
  int<lower=1> J; // number of legislators
  int<lower=1> K; // number of bills
  int<lower=1> N; // number of observations
  int<lower=1,upper=J> j[N]; // legislator for observation n
  int<lower=1,upper=K> k[N]; // bill for observation n
  int<lower=0,upper=1> y[N]; // vote of observation n
  real party[J]; // party of legislator j (0 for D/I, 1 for R)
}
parameters {
  real alpha[K];             
  real beta[K];   
  real theta[J]; # realized ideology
  real gamma[J]; # unobserved, true ideology
  real beta_party; # effect of party ID
}
model {
  alpha ~ normal(0, 5); 
  beta ~ normal(0, 5);
  beta_party ~ normal(0, 2);
  for (i in 1:J){
    theta[i] ~ normal(gamma[i], 1); # hierarchical structure for obs. ideol.
    gamma[i] ~ normal(beta_party * party[i], 1); # effect of party ID
  };
  for (n in 1:N)
    y[n] ~  bernoulli_logit( theta[j[n]] * beta[k[n]] - alpha[k[n]] );
}
'

repub <- ifelse(rc$legis.data$party=="R", 0, 1)
stan.data <- list(J=J, K=K, N=N, j=j, k=k, y=y, party=repub)

stan.fit <- stan(model_code=stan.code, data=stan.data, iter=500, warmup=200,
                 chains=1, thin=2)

monitor(stan.fit)


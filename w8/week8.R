################################################################
## Panel Data
## Author: Jason Qiang Guo
## Quant III Week 8
## November 3rd
################################################################
library(foreign)
library(plm)
library(lme4)
library(pglm)
library(survival)

library(MSwM)

setwd("~/Dropbox/Quant3_TA/TA/w8")
##################
## Linear Model ##
##################

# linear model fixed effects in R
data("Produc", package = "plm")
fe <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp,
          data = Produc, index = c("state","year"), model = "within")
summary(fe)

# random effects 
re <- plm(log(gsp) ~  log(pcap) + log(pc) + log(emp) + unemp,
          data = Produc, index = c("state","year"), model = "random")
summary(re)

# between effects 
be <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp,
          data = Produc, index = c("state","year"), model = "between")
summary(be)

# another way to do it: lmer function is better for varying intercept and varying coefficients
re2 <- lmer(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp + (1 | state),
          data = Produc)
summary(re2)

# hausman test, p-value is small means we reject the null hypothesis that random effects is appropriate
phtest(fe, re)


# Incidental Parameter Problem
# Estimation of standard errors are likely to be inconsistent when T is small while N is large
# for probability model dealing with panel data due to the none-linearity nature of the link function

# but still we can do it with conditional logit model

BESpanel <- read.csv(url("http://www.polsci.ucsb.edu/faculty/glasgow/ps206/ps206data6.txt"), header=T, sep="\t")

# conditional logit model with fixed effects
logit.fe <- clogit(natpess ~ unsatdem+attpol+infpol+age + strata(besid), data=BESpanel)
summary(logit.fe)

# using cox model with single discrete time, exactly the same thing!
coxph(Surv(time = rep(1, 16186), natpess) ~ unsatdem + attpol+infpol+age + strata(besid), data=BESpanel, method = "exact")

# markov switching model
y <- rep(NA, 300)
# state 1
for (i in c(1:100, 151:180, 251:300)){y[i] <- -3 * x[i] + 1 + rnorm(1)}
for (i in c(101:150, 181:250)){y[i] <- 2 * x[i] + 3 + rnorm(1)}
plot(ts(y))

mod <- lm(y~x)
mod.msm <- msmFit(mod, k = 2, sw = c(T, T, F))
plotProb(mod.msm, which=1)

##################################################
## Binary Choice Model 
## Instructor: Jason Guo
## Quant III Lab 2
#################################################

# Probit model
# Actual Data
d2 <- read.csv("week2_data.csv", stringsAsFactors=F)

# probit likelihood function
probit.loglk_manual <- function(beta, X, y){
  loglk <- sum(y * pnorm(X%*%beta, log.p=TRUE) + 
                 (1-y) * pnorm(-(X%*%beta), log.p=TRUE))
  return(-loglk)
}

X <- model.matrix(~age + female + high.school + college + nytimes, d2)
y <- d2$twitter

init <- rnorm(n = ncol(X))
fit <- optim(par = init,fn = probit.loglk_manual, method = "BFGS", X=X, y=y, hessian = T)
se.fit <- sqrt(diag(solve(fit$hessian)))

d3 <- as.matrix(data.frame(c=1, age = 18:80, female=median(d2$female),
                            high.school=median(d2$high.school), 
                            college=median(d2$college),
                            nytimes = median(d2$nytimes)))

# predicted probabilities on age
predicted_p <- pnorm(d3%*%fit$par)

# s.e. of predicted probabilities
p.se <- sqrt(diag(d3 %*% solve(fit$hessian) %*% t(d3)))

# 95% confidence interval
CIupper <- pnorm(d3%*%fit$par + p.se*1.96)
CIlower <- pnorm(d3%*%fit$par - p.se*1.96)

dnorm(d3%*%fit$par) * fit$par[2]
plot(18:80, dnorm(d3%*%fit$par) * fit$par[2], type = "l")

plot(18:80, predicted_p, type = "l", xlab = age, ylab="Pr(Twitter=1)", ylim = c(0, 0.40))
lines(18:80, CIupper, col="grey80")
lines(18:80, CIlower, col="grey80")

# easy way to do it with built-in functions in R
probit <- glm(twitter ~ age + female + high.school + college + nytimes,
              data=d2, family=binomial(link="probit"))

db <- data.frame(age = 18:80, female=median(db$female),
                 high.school=median(db$high.school), 
                 college=median(db$college),
                 nytimes = median(db$nytimes))

pred <- predict(probit, db, type="response", se.fit=TRUE)
plot(18:80, pred$fit, type="l", xlab="age", ylab="Pr(Twitter=1)", ylim=c(0, .40))
lines(18:80, pred$fit - pred$se.fit * 1.96, col="grey80")
lines(18:80, pred$fit + pred$se.fit * 1.96, col="grey80")

plot.data <- data.frame(age=18:80, pred = pred$fit, lo = pred$fit - pred$se.fit * 1.96, hi=pred$fit + pred$se.fit * 1.96)

# much fancier plot with ggplot2
library(ggplot2)
plot <- ggplot(plot.data, aes(x=age, y=pred)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.1, fill="red") + 
  geom_line() +
  scale_x_continuous(limits=c(20,80)) +
  scale_y_continuous("Predicted Prob. of Twitter Use", limits=c(0, .35)) +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) 
plot

# marginal effects plot
xb <- d3 %*% fit$par
plot(18:80, dnorm(xb) * fit$par[2], type = "l", xlab = "age", ylab = "Marginal Effects of Age on Probability of Twitter Use")
lines(18:80, dnorm(xb) * fit$par[2] + 1.96 * dnorm(xb)* se.fit[2], col = "grey80")
lines(18:80, dnorm(xb) * fit$par[2] - 1.96 * dnorm(xb)* se.fit[2], col = "grey80")

## Identification of Binary Choice Model
## What happens when we try to estimate a logit/probit model with intercept
## and a threshold parameter?

## log.lk function adding an ancillary parameter
logit.loglk <- function(pars, X, y){
  a <- pars[1]
  b <- pars[2:length(pars)]
  loglk <- sum(y * plogis(a + X%*%b, log.p=TRUE) + (1-y) * plogis(-(a + X%*%b), log.p=TRUE))
  return(-loglk)
}

# data
y <- d2$twitter
X <- model.matrix(~age+female, d2)
k <- dim(X)[[2]] ## number of variables



# 'correct' values
summary(glm(twitter ~ age + female,
            data=d2, family=binomial(link="logit")))

# ML estimation
initial <- c(0, 0, 1, 1)
logit <- optim(par = initial, logit.loglk, X=X, y=y, method="BFGS", hessian = T)
logit$par

# what happens when we choose different initial values?
initial <- c(-10, 10, rep(0, k-1))
logit <- optim(par = initial, logit.loglk, X=X, y=y, method="BFGS", hessian = T)
logit$par

# what happens to standard errors?
initial <- c(0, 0, 1, 1)
logit <- optim(par = initial, logit.loglk, X=X, y=y, method="BFGS", hessian = T)
logit$par
round(sqrt(diag(solve(logit$hessian))), 3)

initial <- c(-10, 10, rep(0, k-1))
logit <- optim(par = initial, logit.loglk, X=X, y=y, method="BFGS", hessian = T)
logit$par
round(sqrt(diag(solve(logit$hessian))), 3)

intercepts <- seq(-10, 10, .10)
thresholds <- seq(-10, 10, .10)
lks <- matrix(NA, nrow=length(intercepts), ncol=length(thresholds))


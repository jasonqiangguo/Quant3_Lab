library(foreign)
library(survival)
install.packages("KMsurv")
library(KMsurv)
library(Zelig)

#Independent variables:
  
#durat: Duration of the government
#invest: Investiture dummy variable - existence of a legal requirement for legislative investiture and is a hurdle that should diminish average duration by causing some governments to fail very quickly
#polar: Polarization index - measure of support for extremist parties, indicating bargaining system complexity and diminished duration. varies 0-43, with a mean of 15.3
#fract: Fractionalization - index that characterizes the number and size of parties in parliment; increased complexity is expected to reduce cabinet duration. varies 349-868, with a mean of 719.
#numst: Numerical status - dummy variable that distinguishes between majority (coded 1) and minority (coded 0) governments, with majority governments expected to last longer
#format: Formation attempts - the number of attempts to form a government during the crisis. the more foiled attempts, the more complex the bargaining environment, and the shorter the cabinet duration. varies 1-8, with a mean of 1.9
#postelec: Postelection - modeled as a dummy variable coded 1 if the government formed immediately after the election, and 0 otherwise. forming in midterm may indicate situational instability not picked up by other variables
#caretakr: Caretaker government - control for caretaker governments of shorter durations that hold office while more 'permanent' replacement is negotiated
#failure: a cabinet was dissolved,1; censored,0

cab <- read.dta("cabinet.dta")

# Nonparametric estimation, Kaplan-Meier estimator
nonpar <- survfit(Surv(durat, event= failure) ~ 1, data=cab)
plot(nonpar, xlab="Months", ylab="F(t)", fun="event")

# (note that first we need to fit a model with just a constant)

# survivor function S(t) = 1 - F(t) OR CDF for (1-f(t))
plot(nonpar, xlab="Months", ylab="S(t)")

# hazard rate h(t) = f(t) / S(t)  OR  f(t) / (1 - F(t))
# discrete hazard rate (note that Stata applies some smoothing to this)
library(muhaz)
fit2 <- kphaz.fit(cab$durat, cab$failure)
kphaz.plot(fit2)

#############     PARAMETRIC MODELS      #######################

################################################################
# EXPONENTIAL DISTRIBUTION: CONSTANT HAZARD RATE
################################################################

# Estimation
exponential <- survreg(Surv(durat, failure) ~ 1 + invest + polar + fract + numst + format + 
  caretakr, dist = "exponential", data = cab)

summary(exponential)
round(exponential$coefficients, 3)

## Use Zelig to Simulate Expected Duration Time E(T) for polar Profile ##
exponential.zelig <- zelig(Surv(durat, failure) ~ invest + polar + fract + numst + format + 
        caretakr, model = "exp", data = cab, cite = F)
polar.seq <- seq(min(cab$polar), max(cab$polar), length = 50)
x.polar <- setx(exponential.zelig, invest = 0, polar = polar.seq, fract = mean(cab$fract), 
                numst = 1, format = mean(cab$format), caretakr = 0)
sim.polar <- sim(exponential.zelig, x = x.polar)

pe.polar <- rep(0, 50)
for (i in 1:50){
  pe.polar[i] <- apply(sim.polar$getqi(qi="ev", xvalue="range")[[i]], 2, mean)
}

# 95% Confidence interval
lo.polar <- rep(0, 50)
for (i in 1:50){
  lo.polar[i] <- apply(sim.polar$getqi(qi="ev", xvalue="range")[[i]], 2, quantile, prob = 0.025)
}

hi.polar <- rep(0, 50)
for (i in 1:50){
  hi.polar[i] <- apply(sim.polar$getqi(qi="ev", xvalue="range")[[i]], 2, quantile, prob = 0.975)
}

# Make the plot
par(mar = c(4, 6, 0.1, 0.5))
plot(polar.seq, pe.polar, type = "n", xlab = "", ylab = "", ylim = c(0, 70), 
     axes = FALSE)
abline(v = seq(min(polar.seq), max(polar.seq), length = 10), col = "gray75", 
       lty = 3)
abline(h = seq(0, 70, by = 5), col = "gray75", lty = 3)
lines(polar.seq, pe.polar, lwd = 3, lty = 1)
lines(polar.seq, lo.polar, lwd = 2, lty = 2)
lines(polar.seq, hi.polar, lwd = 2, lty = 2)
title(ylab = expression("Expected Cabinet Duration"), line = 4, cex.lab = 1.5)
title(xlab = expression("Support for Extremist Parties"), line = 2.75, cex.lab = 1.5)
axis(1)
axis(2, at = seq(0, 70, by = 5), las = 2, cex.axis = 1.1)
box()
rug(jitter(cab$polar), ticksize = 0.015)
legend("topright", bty = "n", c(expression("Point Estimate"), expression("95% Conf. Interval")), 
       lty = c(1, 2), lwd = c(3, 2), cex = 1.25)

# Marginal effects
x <- c(1, median(cab$invest), mean(cab$polar), mean(cab$fract), 1, mean(cab$format), median(cab$caretakr))
b <- exponential$coefficients
b[5] * exp(x %*% b)

# baseline hazard rate
t <- seq(0, 70, 1)
lambda.i <- exp(-predict(exponential, cab, type="linear"))
lambda <- mean(lambda.i, na.rm=TRUE)

hazard <- lambda
plot(t, rep(hazard, length(t)), type="l", main="Exponential", xlab="Months",
     ylab="Hazard Rate")


################################################################
# WEIBULL DISTRIBUTION: MONOTONICALLY INCREASING/DECREASING HAZARD RATE
################################################################

weibull <- survreg(Surv(durat, failure) ~ 1 + invest + polar + fract + numst + format + 
                     caretakr, dist = "weibull", data = cab)

summary(weibull)
round(weibull$coefficients, 3)

## Use Zelig to Simulate Expected Duration Time E(T) for polar Profile ##
weibull.zelig <- zelig(Surv(durat, failure) ~ invest + polar + fract + numst + format + 
                             caretakr, model = "weibull", data = cab, cite = F)
polar.seq <- seq(min(cab$polar), max(cab$polar), length = 50)
x.polar <- setx(weibull.zelig, invest = 0, polar = polar.seq, fract = mean(cab$fract), 
                numst = 1, format = mean(cab$format), caretakr = 0)
sim.polar <- sim(weibull.zelig, x = x.polar)

pe.polar <- rep(0, 50)
for (i in 1:50){
  pe.polar[i] <- apply(sim.polar$getqi(qi="ev", xvalue="range")[[i]], 2, mean)
}

# 95% Confidence interval
lo.polar <- rep(0, 50)
for (i in 1:50){
  lo.polar[i] <- apply(sim.polar$getqi(qi="ev", xvalue="range")[[i]], 2, quantile, prob = 0.025)
}

hi.polar <- rep(0, 50)
for (i in 1:50){
  hi.polar[i] <- apply(sim.polar$getqi(qi="ev", xvalue="range")[[i]], 2, quantile, prob = 0.975)
}

# Make the plot
par(mar = c(4, 6, 0.1, 0.5))
plot(polar.seq, pe.polar, type = "n", xlab = "", ylab = "", ylim = c(0, 70), 
     axes = FALSE)
abline(v = seq(min(polar.seq), max(polar.seq), length = 10), col = "gray75", 
       lty = 3)
abline(h = seq(0, 70, by = 5), col = "gray75", lty = 3)
lines(polar.seq, pe.polar, lwd = 3, lty = 1)
lines(polar.seq, lo.polar, lwd = 2, lty = 2)
lines(polar.seq, hi.polar, lwd = 2, lty = 2)
title(ylab = expression("Expected Cabinet Duration"), line = 4, cex.lab = 1.5)
title(xlab = expression("Support for Extremist Parties"), line = 2.75, cex.lab = 1.5)
axis(1)
axis(2, at = seq(0, 70, by = 5), las = 2, cex.axis = 1.1)
box()
rug(jitter(cab$polar), ticksize = 0.015)
legend("topright", bty = "n", c(expression("Point Estimate"), expression("95% Conf. Interval")), 
       lty = c(1, 2), lwd = c(3, 2), cex = 1.25)

# Marginal effects
x <- c(1, median(cab$invest), mean(cab$polar), mean(cab$fract), 1, mean(cab$format), median(cab$caretakr))
b <- weibull$coefficients
b[5] * exp(x %*% b)

# Hazard rate
lambda.i <- exp(-predict(weibull, cab, type="linear"))
lambda <- mean(lambda.i, na.rm=TRUE)
t <- seq(0,70,1)
p <- 1/weibull$scale
scale <- weibull$scale
hazard <- lambda * p * (lambda * t)^(p-1)
plot(t, hazard, type="l", main="Weibull", 
     xlab="Months", ylab="Hazard Rate")


################################################################
# LOGNORMAL DISTRIBUTION: NONMONOTONIC  HAZARD RATE
################################################################

lognormal <- survreg(Surv(durat, failure) ~ 1 + invest + polar + fract + numst + format + 
                       caretakr, dist = "lognormal", data = cab)

summary(lognormal)
round(lognormal$coefficients, 3)

# Marginal effects
x <- c(1, median(cab$invest), mean(cab$polar), mean(cab$fract), 1, mean(cab$format), median(cab$caretakr))
b <- lognormal$coefficients
b[5] * exp(x %*% b)

# Hazard Rate in Lognormal
lambda.i <- exp(-predict(lognormal, cab, type="linear"))
lambda <- mean(lambda.i, na.rm=TRUE)
p <- 1/exponential$scale
pdf <- (2*pi)^{-1/2} * p * t^{-1} * exp((-p^2 * (log(lambda*t))^2)/2)
cdf <- 1 - pnorm(p*log(lambda*t))

hazard <- pdf/cdf
plot(t, hazard, type="l", main="Log-normal", 
     xlab="Months", ylab="Hazard Rate")



################################################################
# COX MODEL
################################################################

data2 <- read.dta("civil_cox.dta")
names(data2)[36:37] <- c("t", "t0")

(cox.model <- coxph(Surv(date0, date1, event=cens, type="counting") ~ 
                      gini_m + ginmis + rgdpch + elf + elf2 + logpop + y70stv + y80stv + 
                      y90stv + d2 + d3 + d4 + cluster(indsp), data=data2, method="efron"))

hazard <- basehaz(cox.model)
plot(hazard$time, hazard$hazard, type="l", ylab="baseline hazard", xlab="time")

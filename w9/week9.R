################################################################
## Selection Model
## Author: Jason Qiang Guo
## Quant III Week 9
## November 10th
################################################################

pkg <- c("foreign","VGAM", "dplyr", "sampleSelection")

lapply(pkg, require, character.only = TRUE)

data("Mroz87")
Mroz87$kids <- (Mroz87$kids5 + Mroz87$kids618 > 0)
head(Mroz87)

# Female labor supply (lfp = labour force participation)

## Outcome equations without correcting for selection

## Comparison of linear regression and selection model
# ols regression 
outcome1 <- lm(wage ~ exper, data = Mroz87)
summary(outcome1)

# heckman two-step estimation
selection1 <- selection(selection = lfp ~ age + I(age^2) + faminc + kids + educ, outcome = wage ~ exper, 
                        data = Mroz87, method = "2step")
summary(selection1)

plot(Mroz87$wage ~ Mroz87$exper)
curve(outcome1$coeff[1] + outcome1$coeff[2]*x, col="black", lwd="2", add=TRUE)
curve(selection1$coeff[1] + selection1$coeff[2]*x, col="orange", lwd="2", add=TRUE)


## A more complete model comparison
# ols regression
outcome2 <- lm(wage ~ exper + I( exper^2 ) + educ + city, data = Mroz87)
summary(outcome1)

## Correcting for selection
# two-step estimation
selection.twostep2 <- selection(selection = lfp ~ age + I(age^2) + faminc + kids + educ, outcome = wage ~ exper + I(exper^2) + educ + city, 
                                data = Mroz87, method = "2step")
summary(selection.twostep2)

# full likelihood estimation
selection.mle <- selection(selection = lfp ~ age + I(age^2) + faminc + kids + educ, outcome = wage ~ exper + I(exper^2) + educ + city, 
                           data = Mroz87, method = "mle")
summary(selection.mle)


# Heckman model selection "by hand" #

seleqn1 <- glm(lfp ~ age + I(age^2) + faminc + kids + educ, family=binomial(link="probit"), data=Mroz87)
summary(seleqn1)

# Calculate inverse Mills ratio by hand #

Mroz87$IMR <- dnorm(seleqn1$linear.predictors)/pnorm(seleqn1$linear.predictors)

# Outcome equation correcting for selection #

outeqn1 <- lm(wage ~ exper + I(exper^2) + educ + city + IMR, data=Mroz87, subset=(lfp==1))
summary(outeqn1)

# compare to selection function using MLE -- coefficients right, se's wrong
summary(selection.mle)


# interpretation
# If our independent variables does not appear in the selection equation, we can interpret beta as in linear regression
# If it does appear in the selection equation, we must calculate:

beta.educ.sel <- selection.twostep2$coefficients[6]
beta.educ.out <- selection.twostep2$coefficients[10]
beta.IMR <- selection.twostep2$coefficients[12]
delta <- selection.twostep2$imrDelta

marginal.effect <- beta.educ.out - beta.educ.sel * beta.IMR * delta 
mr2 <- marginal.effect * Mroz87$educ


plot(Mroz87$wage ~ Mroz87$educ)
lines(mr2 ~ Mroz87$educ, type="l", col="green", lwd="2")

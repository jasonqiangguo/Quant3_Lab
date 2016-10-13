library(foreign)
library(ggplot2)
library(MASS)
library(mlogit)
library(dplyr)

setwd("~/Dropbox/Quant3_TA/TA/w4")

data <- read.dta("dutch271cond.dta")
str(data)

## data cleaning: delete rows of the individuals for whom values of some of their alternative specific variables are missing
for (i in 1:dim(data)[2]) {
  index <- which(data$respid %in% data$respid[is.na(data[,i]) == TRUE], )
  if (length(index) == 0) {
    data <- data    
  }
  else {
    data <- data[-index,]
  }
}

#vote: VOTE indicates for each line 1 for pvda, 2 for cda, 3 for vvd and 4 for d66
TM <- mlogit.data(data, choice = "choice", shape = "long", chid.var = "respid", alt.levels = c("pvda", "cda", "vvd", "d66"))

# policy distance (alternative specific) variables: disabor, disnuke

# pure conditional model
fA <- mFormula(choice ~ disabor + disnuke|0|0)
multiA <- mlogit(fA, TM, reflevel = "pvda")
summary(multiA)

# we dont have a constant term in the estimation, because if a constant is added to the utility of all alternatives, then the alternative with the highest
# utility does not change and the effect of constant term cancels each other out and therefore the constant becomes unidentifiable. Similarly, the level of
# utility does not matter for the analyst. The fact that only differences in utility matter has implications for the identification of discrete choice models.
# In particular, it means that the only parameters that can be estimated are those that captures differences across alternatives. Unless we add 
# alternative-specific constants, we should not see a constant term.

# drop party 4 and individuals that choose party 4
# labor consuming way of doing this
id_v4 <- data$respid[data$choice == 1 & data$vote == 4]
index <- which(data$respid %in% id_v4)
data_dpv4 <- data[-index,]
data_dpv4 <- filter(data_dpv4, vote != 4)

# reshape data
TM_dp_v4 <- mlogit.data(data_dpv4, choice = "choice", shape = "long", chid.var = "respid", alt.levels = c("pvda", "cda", "vvd"))

# policy distance (alternative specific) variables: disabor, disnuke

fA_dp_v4 <- mFormula(choice ~ disabor + disnuke|0|0)
multiA_dp_v4 <- mlogit(fA_dp_v4, TM_dp_v4, reflevel = "pvda")
summary(multiA_dp_v4)

# or, a much easier way to drop alternatives in estimation by restricting the selection of alternatives with alt.subset
multiA_dp_v4_2 <- mlogit(fA, TM, reflevel = "pvda", alt.subset = c("pvda", "cda", "vvd"))
summary(multiA_dp_v4_2)
# The result shows that dropping a choice does not change the estimation, especially if we take a look at the frequencies of alternatives,
# the ratio does not change by a lot

# hausman-macfadden test
# built-in function 
hmftest(multiA, multiA_dp_v4_2)

# or by hand
coeffF1 <- as.matrix(multiA$coefficients)
coeffF <- coeffF1[order(rownames(coeffF1))]
coeffF <- coeffF[1:2]
coeffR1 <- as.matrix(multiA_dp_v4_2$coefficients)
coeffR <- coeffR1[order(rownames(coeffR1))]
coeffD <- coeffR - coeffF

covmatF1 <- - solve(multiA$hessian)
covmatF <- covmatF1[order(rownames(covmatF1)),order(colnames(covmatF1))]
covmatF <- covmatF[1:2,1:2]
covmatR1 <- - solve(multiA_dp_v4_2$hessian)
covmatR <- covmatR1[order(rownames(covmatR1)),order(colnames(covmatR1))]
covmatD <- covmatR - covmatF

HMF <- t(coeffD)%*%solve(covmatD)%*%coeffD
HMF

# results should be the same

# add sociological variables (individual specific variables)
fa_ind <- mFormula(choice ~ disabor + disnuke | age + income)
multiA_ind <- mlogit(fa_ind, TM, reflevel = "pvda")
summary(multiA_ind)

## add sociological variables and also individuals' LR to make the model unconditional
fa_uncond <- mFormula(choice ~ 0 | abor + nuke + age + income)
multiA_uncond <- mlogit(fa_uncond, TM, reflevel = "pvda")
summary(multiA_uncond)

## plot odds ratio
plot.data <- data.frame(exp(cbind(OR = coef(multiA_uncond_cda), confint(multiA_uncond_cda))))[-(1:3),]
plot.data$party <- rep(c("pvda","vvd","d66"), 4)
plot.data$vname <- c(rep("abor", 3), rep("nuke", 3), rep("age", 3), rep("income", 3))
names(plot.data) <- c("OR","lb", "ub", "party", "vname")

ggplot(plot.data, aes(x = party, y = OR, ymin = lb, ymax = ub)) + geom_pointrange(aes(col = factor(vname)), position=position_dodge(width=0.30)) + 
  ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 0)) + scale_color_discrete(name = "vname") + xlab("")


## change the base category
multiA_uncond_cda <- mlogit(fa_uncond, TM, reflevel = "cda")
summary(multiA_uncond_cda)

## what if we have a dataset in which each row is an observation like in 271uncond.dta?
## one way is to directly run a multinomial logit regression using the function multinom()
data_uncond <- read.dta("dutch271uncond.dta")
mnformula <- mFormula(party ~ abor + nuke + inc + left + age + relig)

library(nnet)
mnl <- multinom(mnformula, data_uncond)
mnl

## the other way is to convert data such that function mlogit can execute, remember the multinomial logit
## regression is equivalent to unconditional logit regression as all the variables are individual specific
## and therefore we fit a coefficient for each alternative

# reshape the data and shape = "wide" means the original dataset treats each observation as a row
mlogitSample <- mlogit.data(data_uncond, choice = "party", shape = "wide")

mnlformula <- mFormula(party ~ 0|abor + nuke + inc + left + age + relig)
mnlogit <- mlogit(mnlformula, mlogitSample, reflevel = "0")
summary(mnlogit)

# the results should be the same


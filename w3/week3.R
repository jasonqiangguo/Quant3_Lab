#########################################################
## Count Data Modeling and Multiple Category Data Modeling
## Instructor: Jason Guo
## Quant III Lab 3
#########################################################

library(foreign)
setwd("/Users/qiangguo/Dropbox/Quant3_TA/TA/w3")

# count data with overdispersion
data <- read.dta("soccer_data.dta")

# Miguel, Edward, Sebasti?n M. Saiegh, and Shanker Satyanath. 
# "Civil war exposure and violence." Economics & Politics 23.1 (2011): 59-73.

library(ggplot2)
ggplot(data, aes(yellow_card)) + geom_histogram()
## estimating poisson model
#---------------------------------------------------------------------
summary(
  poisson <- glm(yellow_card ~ civwar + age, family = "poisson", data = data)
)

# calculate overdispersion parameter
sqrt(sum((residuals(poisson, type="pearson"))^2)/5417)


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
     xlab="Age", ylab="Poisson: Marginal effect of age", ylim=c(0.02, 0.07))
lines(age, marg[2,], col="grey")
lines(age, marg[3,], col="grey")

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
     xlab="Age", ylab="NB: Marginal effect of age", ylim=c(0.02, 0.07))
lines(age, marg[2,], col="grey")
lines(age, marg[3,], col="grey")

# but how different is it from what we would expect from the proper dist?
xb=cbind(1,data$civwar,data$age)%*%negbin$coefficients

zeros=numeric(1000)
for (i in 1:1000)
{
  zeros[i]=sum(rpois(n = length(xb),lambda = exp(xb))==0)
}

plot(density(zeros))

plot(density(zeros),xlim = c(400,2000))
abline(v = sum(data$yellow_card==0))

## let us try using zero-inflated model
#---------------------------------------------------------------------
library(pscl)
zip <- zeroinfl(yellow_card ~ civwar + age, data=data)
summary(zip)

zip.coefs<- zip$coefficients
zip.count <- zip.coefs$count
zip.zero <- zip.coefs$zero
zip.se <- sqrt(diag(vcov(zip)))


############################################
# ordered probit model
############################################

# reading the data
data <- read.csv("week3_data.csv", stringsAsFactors=F)


# You can do this manually
# Code to run a four-category ordered probit
# Chris Adolph (University of Washington)

# Likelihood for 4 category ordered probit
llk.oprobit4 <- function(param, x, y) {
  # preliminaries
  os <- rep(1, nrow(x))
  x <- cbind(os, x)  
  b <- param[1:ncol(x)]
  t2 <- param[(ncol(x)+1)]
  t3 <- param[(ncol(x)+2)]
  
  # probabilities and penalty function
  xb <- x%*%b
  
  p1 <- log(pnorm(-xb))
  
  # penalty function to keep t2>0
  if (t2<=0)  p2 <- -(abs(t2)*10000) else p2=log(pnorm(t2-xb)-pnorm(-xb))
  
  # penalty to keep t3>t2
  if (t3<=t2) p3 <- -((t2-t3)*10000) else p3=log(pnorm(t3-xb)-pnorm(t2-xb))     
  
  p4 <- log(1-pnorm(t3-xb)) 
  
  # -1 * log likelihood (optim is a minimizer)
  -sum(cbind(y==1,y==2,y==3,y==4) * cbind(p1,p2,p3,p4))
}

# Preparing data

y <- data$nytimes   # How often respondent reads NY Times
## 4 = regularly
## 3 = sometimes
## 2 = hardly ever
## 1 = never

x <- cbind(data$age, data$female, data$high.school, data$college)

# Use optim directly
ls.result <- lm(y~x)                    # use ls estimates as starting values
stval <- c(ls.result$coefficients,2,3)  # initial guesses
oprobit.result <- optim(stval, llk.oprobit4, method="BFGS", x=x, y=y, hessian=T)

pe <- oprobit.result$par                # point estimates
vc <- solve(oprobit.result$hessian)     # var-cov matrix
se <- sqrt(diag(vc))                    # standard errors

results <- cbind(pe, se)
dimnames(results) <- list(
  c("Intercept", "age", "female", "high.school", "college", 
    "tau2", "tau3"),
  c("Value", "Std. Error"))
round(results, 5)


############################################################################
# or, you can use R built-in function
# Use MASS polr to do ordered probit
############################################################################
library(MASS)
oprobit <- polr(
  factor(nytimes) ~ age + female + high.school + college, 
  data=data, method="probit", Hess = TRUE)
summary(oprobit)
oprobit$Hessian
oprobit$zeta
oprobit$coefficients

# plot predicted probabilities 
newdata <- data.frame(age = 18:80, female=mean(data$female),
                      high.school=mean(data$high.school), college=mean(data$college))
pred <- predict(oprobit, newdata, type="probs", se.fit = TRUE)
str(pred)
plot(18:80, pred[,3], type="l")


# using simulation to get standard errors
set.seed(123)
coefs <- mvrnorm(n=10000, mu=c(oprobit$coefficients, oprobit$zeta), Sigma=solve(oprobit$Hessian))

sim.results <- matrix(data = NA,nrow = 63,ncol = 10000)
age=18:80
for (i in 1:length(age)){
  
  x <- c(age[i], mean(data$female), 
         mean(data$high.school), mean(data$college))
  xb <- x %*% t(coefs[,1:4])
  # recall that pr(y=3)=Phi(tau3-xb)-Phi(tau2-xb)
  tau2 <- coefs[,6]; tau3 <- coefs[,7]
  probs <- pnorm(tau3-xb) - pnorm(tau2-xb)
  sim.results[i,] <- probs
}


# getting quantiles
plot.data <- apply(sim.results, 1, function(x) c(quantile(x, .025), mean(x), quantile(x, .975)))

plot(18:80, plot.data[2,], type="l", ylim=c(.05, .25),
     xlab="age", ylab="Pr(y=3)")
lines(18:80, plot.data[1,], lty=3)
lines(18:80, plot.data[3,], lty=3)


# use built-in function to run multinomial logit model
rm(list=ls())

load(file = "united_fmly_2012.rda")

attach(fmly2012)
# 1 Married
# 2 Widowed
# 3 Divorced
# 4 Separated
# 5 Never married
table(marital1)

# form=as.formula(
#   marital1~ref_race+sex_ref+bls_urbn+region+cutenure+age_ref+educ_ref)

fmly2012$sex_ref=as.numeric(sex_ref)

form=as.formula(marital1~sex_ref+age_ref)

library(nnet)
model=multinom(formula = form,data = fmly2012,Hess = T,maxit=1000)

summary(model)

summed_model <- summary(model)

# par(mfrow=c(1,1))

newdata_1 <- data.frame(sex_ref=1,age_ref=median(age_ref))
p_ref1 <- predict(object <- model,newdata = newdata_1,type="probs")
plot(x = 1:5,y=predict(object = model,newdata = newdata_1,type="probs"),
     type="h",ylab = "pie",xlab="1")

newdata_2 <- data.frame(sex_ref=2,age_ref=median(age_ref))
p_ref2 <- predict(object = model,newdata = newdata_2,type="probs")
plot(x = 1:5,y=predict(object = model,newdata = newdata_2,type="probs"),
     type="h",ylab = "pie",xlab="2")

dim(summed_model$coefficients)
dim(summed_model$standard.errors)
sqrt(diag(solve(model$Hessian)))

H=solve(model$Hessian)
# mat=matrix(data = 1,nrow = dim(H)[1],ncol = dim(H)[2])
# diag(mat)=diag(solve(model$Hessian))



ITER=1000
library(MASS)
sim_mat=mvrnorm(n=ITER, mu=as.numeric(summed_model$coefficients), 
                Sigma=H)

martial=c("1","2","3","4","5")
pie0=matrix(data = NA,ncol = length(martial),nrow = ITER)
pie1=matrix(data = NA,ncol = length(martial),nrow = ITER)

for (k in 1:ITER)
{
  B=matrix(data = sim_mat[k,],nrow = 4)
  
  #   as.numeric(sim_mat[ITER,])
  
  
  #   B_=as.matrix(rbind(0,summed_model$coefficients))
  B=as.matrix(rbind(0,B))
  
  X_0=as.matrix(cbind(1,newdata_1))
  X_1=as.matrix(cbind(1,newdata_2))
  deno_0=sum(exp(X_0%*%t(B)))
  deno_1=sum(exp(X_1%*%t(B)))
  
  for (j in 1:ncol(pie0))
  {
    pie0[k,j]=exp(X_0%*%B[j,])/deno_0
    pie1[k,j]=exp(X_1%*%B[j,])/deno_1
  }
  
}


boxplot(pie0)
boxplot(pie1)


# Now we want to plot the change by an individual's gender in his/her predicted probabilites for each choice 
library(ggplot2)
plot.data <- data.frame(1, 2, p_ref1, p_ref2)
Marital_Status <- c("Married", "Widowed", "Divorced", "Separated", "Never married")
b <- ggplot(plot.data, aes(X1, p_ref1)) +
  geom_point() + geom_point(aes(X2, p_ref2))
b <- b + geom_segment(aes(xend = X2, yend = p_ref2, color = legend_label))
b <- b + scale_x_continuous(name = "Sex Reference", limits = c(0.8, 2.2), breaks = c(1, 2), labels = c("male", "female"))
b <- b + ylab("Predicted Probilities")
b

# there is another option: we can use zelig to help us draw the predicted probabilities
library(Zelig)
data(turnout)
z.out <- zelig(vote ~ race + educate + age + I(age^2) + income,
               model = "logit", data = turnout)
age.range <- 18:95
x.low <- setx(z.out, educate = 12, age = age.range)
x.high <- setx(z.out, educate = 16, age = age.range)
s.out <- sim(z.out, x = x.low, x1 = x.high)

ci.plot(s.out, xlab = "Age in Years",
        ylab = "Predicted Probability of Voting",
        main = "Effect of Education and Age on Voting Behavior",  ci = c(95, 99, 99.9))
legend(45, 0.52, legend = c("College Education (16 years)",
                            "High School Education (12 years)"), col = c("blue","red"), 
       lty = c("solid"))

#####################################################################
## conditional logit model: the estimation of the multinomial logit
## models with individual and alternative specific variables
#####################################################################

library(mlogit)
# Let us use TravelMode datast as an example. This dataset has both individual specific variables (income and homesize)
# and alternative specific variables like monetary cost (vcost) and travel time (travel). In this case,
# the alternative specific variable vcost is associated with a generic coefficient, but the other alternative specific variable
# travel is associated with an alternative specific coefficient as travel time, although heterogenous across individuals, mostly
# is largely determined by the mode of transportation rather than individuals. As for individual specific variables income and size, 
# we use alternative specific coefficinets for estimation, since given an individual's income and homesize his satisfaction (utility)
# is also determined by the mode of transportation.

data("TravelMode", package = "AER")
head(TravelMode)


# data preparation: here shape = "long" means each row is an alternative.
# The name of the variable that contains the information about the choice situations can
# be indicated using the chid.var argument. alt.var which indicates the name of the
# variable that contains the alternatives
TM <- mlogit.data(TravelMode, choice = "choice", shape = "long", chid.var = "individual", alt.levels = c("air", "train", "bus", "car"))

# There three different types of variables that mlogit can specify
# 1. alternative specific variables x_ij with a generic coefficient \beta,
# 2. individual specific variables z_i with alternative specific coefficient \gamma_j
# 3. alternative specific variables w_ij with an alternative specific coefficient \delta_j .

# you can check the data formula for estimation, and it displays dummy variable-like tricks
# in this formula, three types of variables are separated by "|"
f <- mFormula(choice ~ vcost | income + size | travel)
head(model.matrix(f, TM))

# model estimation
multi <- mlogit(choice ~ vcost | income + size | travel, TM)

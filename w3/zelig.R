par(mfrow=c(1,1))
library(Zelig)
library(MASS)
data(turnout)
z.out <- zelig(vote ~ race + educate + age + I(age^2) + income,
               model = "logit", data = turnout)
age.range <- 18:95
x.low <- setx(z.out, educate = 12, age = age.range)
x.high <- setx(z.out, educate = 16, age = age.range)
s.out <- sim(z.out, x = x.low, x1 = x.high)

s.outout <- sim(z.out, x.low)

ci.plot(s.out, xlab = "Age in Years",
        ylab = "Predicted Probability of Voting",
        main = "Effect of Education and Age on Voting Behavior",  ci = c(95, 99, 99.9))
legend(45, 0.52, legend = c("College Education (16 years)",
                            "High School Education (12 years)"), col = c("blue","red"), 
       lty = c("solid"))
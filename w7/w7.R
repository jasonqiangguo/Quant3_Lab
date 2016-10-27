################################################################
## PCSE
## Author: Jason Qiang Guo
## Quant III Week 7
## October 27th
################################################################

library(foreign)
install.packages("pcse")
library(pcse)
library(car)
library(lmtest)
library(urca)

# load Cross-National panel data on the OECD countries containing information about economic performance,
# government partisanship and labor organization. AGL estimated a model relating political and labor organization variables
# (and some economic controls) to economic growth, unemployment, and inflation. The argument was that economic performance
# in advanced industrialized societies was superior when labor was both encompassing and had political power
# or when labor was weak both in politics and the market. Here we will only look at their model of economic growth.

data(agl)
summary(agl)

# specify the model to test the theory
agl.lm <- lm(growth ~ lagg1 + opengdp + openex + openimp + central +
               +  leftc + inter + as.factor(year), data = agl)

summary(agl.lm)

# d-w test for serial correlation
durbinWatsonTest(agl.lm, max.lag=1)
# b-g test for serial correlation
bgtest(agl.lm, order = 1)
# does not have strong evidence of serial correlation

# Is growth I(1)?
summary(ur.df(agl$growth, type="none", lags=0))

# Is opengpd I(1)?
summary(ur.df(agl$opengdp, type="none", lags=0))

# ...
# do not have strong evidence that any time series is an integrated process

# apply pcse
agl.pcse <- pcse(agl.lm, groupN = agl$country, groupT = agl$year)
summary(agl.pcse)

# We note that the standard error on central has increased a bit, but the standard errors of the other two variables of interest,
# leftc and inter have actually decreased.

####################
# unbalanced panel #
####################

# What if the panel is unbalanced, i.e., some units have shorter time periods than others due to case omission?
# casewise approach: calculate the vcov using the largest balanced subset of the data. 
data(aglUn)
summary(aglUn$country)

aglUn.lm <- lm(growth ~ lagg1 + opengdp + openex + openimp + central +
               +  leftc + inter + as.factor(year), data = aglUn)

summary(aglUn.lm)

agl.pcse <- pcse(aglUn.lm, groupN = aglUn$country, groupT = aglUn$year, pairwise = FALSE)
summary(agl.pcse)

# This warning means that the largest balanced panel only has seven time points, whereas the data runs for 14.

# pairwise approach: calculate the vcov for each pair of units we determine with temporal observations
# overlap between the two.

agl.pcse2 <- pcse(aglUn.lm, groupN = aglUn$country, groupT = aglUn$year, pairwise = TRUE)
summary(agl.pcse2)

# In this case, it is not clear that PCSEs will be correctly estimated although in this case they are 
# not that different from the casewise estimate.


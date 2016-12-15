######################################################
# unsurpervised learning: lasso and ridge regression #
# Jason Guo, Quant 3 Lab                             #
######################################################

pkg <- c("foreign", "ggplot2", "glmnet", "ISLR")
lapply(pkg, require, character.only = TRUE)

set.seed(999)

# load data
Hitters=na.omit(Hitters)
x <- model.matrix(Salary~.-1,data=Hitters) 
y <- Hitters$Salary

# Fitting the model (Lasso: Alpha = 1)
lasso <- glmnet(x, y, alpha = 1)
plot(lasso, xvar = "lambda", label = TRUE)

# Model Selection (cross-validation)
cv.lasso <- cv.glmnet(x, y, alpha=1, type.measure='mse', nfolds = 10)
plot(cv.lasso) # the number on the top of the graph means the number of non-zero coefficients as lambda (penalty) increases
coef(cv.lasso) # non-zero coefficients

# Fitting the model (Ridge: Alpha = 0)
ridge <- glmnet(x, y, alpha = 0)
plot(ridge, xvar = "lambda", label = TRUE)

# Model Selection (cross-validation)
cv.ridge <- cv.glmnet(x, y, alpha=0, type.measure='mse', nfolds = 10)
plot(cv.ridge)
coef(cv.ridge)

######################################################
# unsurpervised learning: CART                       #
# Jason Guo, Quant 3 Lab                             #
######################################################

require(caret)
data <- read.csv("lab3_data.csv")

set.seed(2016)

# create traning set and test set
index_train=createDataPartition(y = data$nytimes, p = .5,list=F)
train=data[index_train,]
test=data[-index_train,]

# set paramesters for cross-validation
control <- trainControl(method="repeatedcv", number=10, repeats=3)

# CART
fit.rpart <- train(factor(nytimes) ~., data=train, method="rpart", metric="Accuracy", trControl=control, tuneLength = 10)
plot(fit.rpart$finalModel)
text(fit.rpart$finalModel)

# Bagged CART
fit.treebag <- train(factor(nytimes) ~., data=train, method="treebag", metric="Accuracy", trControl=control, tuneLength = 10)
predict(fit.treebag, test)

# Random Forest
fit.rf <- train(factor(nytimes)~., data=train, method="rf", metric="Accuracy", trControl=control, tuneLength = 10)
plot(fit.rf)

# summarize results
bagging_results <- resamples(list(treebag=fit.treebag, rf=fit.rf))
summary(bagging_results)
dotplot(bagging_results)



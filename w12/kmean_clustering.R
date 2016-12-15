#############################################
# unsurpervised learning: k-mean clustering #
# Jason Guo, Quant 3 Lab                    #
#############################################

pkg <- c("foreign", "ggplot2", "datasets", "dplyr")
lapply(pkg, require, character.only = TRUE)
head(iris)

d <- dplyr::select(iris, Sepal.Length, Sepal.Width) 

# before clustering
ggplot(d, aes(x = Sepal.Length, y = Sepal.Width)) + geom_point()

# clustering
kc <- kmeans(d,20)

d$kc <- kc$cluster

# after clustering
ggplot(d, aes(x = Sepal.Length, y = Sepal.Width, color = factor(kc))) + geom_point()

# how many clusters do we need?
n.k <- 2:20
average.sum.square <- function(x){
  kc <- kmeans(d, x)
  total <- kc$tot.withinss
  average <- total/x

  return(average)
}

aver.sum.square <- unlist(lapply(2:20, average.sum.square))
plot(aver.sum.square, xlab = "Number of Clusters")
lines(aver.sum.square)


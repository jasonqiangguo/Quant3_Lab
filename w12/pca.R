########################################################
# unsurpervised learning: principal component analysis #
# Jason Guo, Quant 3 Lab                               #
########################################################


pkg <- c("foreign", "ggplot2", "datasets", "dplyr","ggbiplot")
lapply(pkg, require, character.only = TRUE)
head(iris)

d <- dplyr::select(iris, Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) 

# original plot

ir.pca <- prcomp(d, center = TRUE, scale. = TRUE)

# check which variable is most important in each component
sweep(abs(ir.pca$rotation),2,colSums(abs(ir.pca$rotation)),`/`)

# check how many components it takes to explain most of the variance
plot(ir.pca, type = "l")

# plot the first and second PCs 
g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1)
g<- g + theme(legend.direction = 'horizontal',
              legend.position = 'top')
print(g)

mean(ir.pca$x[,4])

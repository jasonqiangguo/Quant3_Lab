#############################################
# unsurpervised learning: mixture models    #
# Jason Guo, Quant 3 Lab                    #
#############################################

pkg <- c("mixtools", "foreign", "ggplot2")
lapply(pkg, require, character.only = TRUE)

set.seed(2016)

N <- 10000
components <- sample(1:2,prob=c(0.3,0.7),size=N,replace=TRUE)
mus <- c(10,7)
sds <- c(1,1)
y <- rnorm(n=N,mean=mus[components],sd=sds[components])

plot(y)
plot(density(y))

# first, lets do it mannually, only consider the mus at this moment and take variance estimate as given
mu_1 <- 0
mu_2 <- 1

pi_1 <- 0.5
pi_2 <- 0.5

for( i in 1:100 ) {
  
  ## Given the observed data, as well as the distribution parameters,
  ## what are the latent variables?
  
  T_1 <- pi_1 * dnorm(y, mu_1)
  T_2 <- pi_2 * dnorm(y, mu_2)
  
  P_1 <- T_1 / (T_1 + T_2)
  P_2 <- T_2 / (T_1 + T_2) ## note: P_2 = 1 - P_1
  
  pi_1 <- mean(P_1)
  pi_2 <- mean(P_2)
  
  ## Given the observed data, as well as the latent variables,
  ## what are the population parameters?
  
  mu_1 <- sum(P_1 * y) / sum(P_1)
  mu_2 <- sum(P_2 * y) / sum(P_2)
  
  ## print the current estimates
  
  print(c(mu_1, mu_2, mean(P_1)))
  
}

# use normalmixEM to fit the model
mix.model <- normalmixEM(y, lambda = 0.3, mu = c(8, 7), sigma = c(1, 1))
mix.model$mu
mix.model$lambda
mix.model$sigma

pred.components <- sample(1:2,prob=mix.model$lambda,size=N,replace=TRUE)
pred.y <- rnorm(n=N,mean=mix.model$mu[components],sd=mix.model$sigma[components])

stack.data <- stack(data.frame(y, pred.y))
ggplot(stack.data, aes(values, color = ind)) + geom_density() +
  scale_color_discrete(name = "group") + theme(panel.grid = element_blank())

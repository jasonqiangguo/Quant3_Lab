#########################################################
## Bayesian Statistics: introduction
## Instructor: Jason Guo
## Quant III Lab 10
#########################################################

################################
# Simulate beta distribution ###
################################
library(foreign)
library(ggplot2)
library(gridExtra)

setwd("~/Dropbox/Quant3_TA/TA/w10")

x <- seq(0, 1, length = 21)
dbeta(x, 1, 1)
pbeta(x, 1, 1)

## Visualization, including limit cases:

pl.beta <- function(a,b, ylim = c(0, 3), xlim = c(0, 1)) {
  x <- seq(0, 1, length = 1025)
  fx <- cbind(dbeta(x, a,b), pbeta(x, a,b))
  f <- fx; f[fx == Inf] <- 1e100
  matplot(x, f, ylab="", type="l", ylim=ylim, xlim = xlim,
          main = sprintf("beta(x, a=%g, b=%g)", a,b))
  abline(h = 0:1, col="gray", lty=3)
  legend("top", c("pdf", "cdf"), col = 1:2, lty = 1:3, bty = "n")
}
pl.beta(1,1)
pl.beta(3,1)
pl.beta(1,3)
pl.beta(0.5,0.5)
pl.beta(2, 4)
pl.beta(2, 2)
pl.beta(3, 7)
pl.beta(0, 0)   ## point masses at  {0, 1}

pl.beta(0, 2)   ## point mass at 0 ; the same as
pl.beta(1, Inf)

pl.beta(Inf, 2) ## point mass at 1 ; the same as
pl.beta(3, 0)

pl.beta(Inf, Inf)# point mass at 1/2



##########################################################################
# beta prior and binomial likelihood, solve for posterior analytically ###
##########################################################################

# Simulate data
n <- 82
s <- 73
y <- sample(rep(c(0, 1), c(n-s, s)))

# since the y follows binomial distribution, \hat{\pi} should be n_successes/n = 73/82 = 0.8902
# let't check with MLE
fn <- function(pi){
  l <- -sum(dbinom(y, 1, pi, log = TRUE))
}

mle <- optim(0.05, fn, method="BFGS", hessian = TRUE)
mle$par

# if prior is an uninformative prior (beta(1,1)), then the posterior should be beta(74, 10)
plot.bayes <- function(a1, b1){
  x <- seq(0, 1, length = 100)
  prior.d <- dbeta(x, a1, b1)
  post.d <- dbeta(x, a1 + s, b1 + n - s)
  data <- data.frame(x, prior.d, post.d)
  p <- ggplot(data, mapping = aes(x = x))  + geom_line(aes(y = post.d, color = "posterior", linetype = "posterior"))
  p <- p + geom_line(mapping = aes(y = prior.d, color = "prior", linetype = "prior"))
  p <- p + theme_bw() + theme(panel.border = element_rect(colour = "black", size = 1),
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              axis.title.y = element_blank(), axis.ticks.y = element_blank(),
                              axis.text.y=element_blank(), legend.key = element_blank())
  p <- p + geom_vline(aes(xintercept = s/n, color = "mle", linetype = "mle"))
  p <- p + xlab(expression(pi)) 
  p <- p + scale_linetype_manual(name = "Lines", values = c("posterior" = 2, "prior" = 2, "mle" = 1))
  p <- p + scale_colour_manual(name = "Lines", values = c("posterior" = "blue", "prior" = "red", "mle" = "black"))
  p <- p + ggtitle(sprintf("beta(%g, %g) prior with n = %g", a1, b1, n))
  p
}

pdf("binomial_beta_2.pdf", height = 8, width = 10)
grid.arrange(plot.bayes(1, 1), plot.bayes(2, 4),
             plot.bayes(2, 12), ncol=2)
dev.off()


#########################################################
## posterior is getting close to MLE when n is large ####
#########################################################

n <- 82 * 20
s <- 73 * 20
y <- sample(rep(c(0, 1), c((n-s)*20, s*20)))

pdf("post_n_large_2.pdf", height = 8, width = 10)
grid.arrange(plot.bayes(1, 1), plot.bayes(2, 4),
             plot.bayes(2, 12), ncol=2)
dev.off()



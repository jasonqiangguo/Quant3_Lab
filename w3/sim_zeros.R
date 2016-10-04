
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


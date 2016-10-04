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

summed_model=summary(model)

par(mfrow=c(1,2))
# par(mfrow=c(1,1))
newdata_1=data.frame(sex_ref=1,age_ref=median(age_ref))
predict(object = model,newdata = newdata_1,type="probs")
plot(x = 1:5,y=predict(object = model,newdata = newdata_1,type="probs"),
     type="h",ylab = "pie",xlab="1")

newdata_2=data.frame(sex_ref=2,age_ref=median(age_ref))
predict(object = model,newdata = newdata_2,type="probs")
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



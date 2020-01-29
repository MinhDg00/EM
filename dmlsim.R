rm(list=ls())
gc(TRUE)
source("JointSol_rep.R")
source("cv_dml.R")

# Matrix Power
"%^%" =  function(x, n) 
  with(eigen(x), vectors %*% (values^n * t(vectors)))

# function to intialize y 
mvsim <- function(x,prob,beta,covar){
  n = dim(x)[1]
  m = ncol(covar)
  y = matrix(0,n,m)
  for (i in 1:n){
    y[i,] = rmvnorm(1,mean=x[i,]%*%beta,sigma=covar)
  }
  y
}
set.seed(13)

n <- 500   # nr of observation 
p <- 80  # nbr of component of x
m  <- 2 # nbr of component of y 
k <- 1  # nbr of mixture model
x <- matrix(rnorm(n*p),n,p) 
prob = 1
beta <- array(0, dim = c(p,m)) 
beta <- rbind(c(1,3),c(3,1),c(1,3),c(1,1),matrix(0, p-4, m))

covar <- array(0, dim = c(m,m))
covar <- matrix(c(2,1,1,2),2,2)


y <- mvsim(x=x,covar=covar,beta=beta,prob=prob)

# Initial for testing line by line
lambda =0.05
Y = y
X = x

#fit dml
fit <- JointSol(y,x,lambda=0.05)
P.est <- fit$P
Phi.est <- fit$Phi

fit
# fit dmlpath
la <- seq(0.5,1.5,length=21)
fdml <- dmlpath(y,x,lambda = la)


apply(fdml$coef, 3, function(x) x[1:10,1])
apply(fdml$coef, 3, function(x) x[1:10,2])

apply(fdml$Phi, 3, function(x) x[1:10,1])
apply(fdml$Phi, 3, function(x) x[1:10,2])

#apply(fdml$P, 3, function(x) x[1:10,1])
#apply(fdml$P, 3, function(x) x[1:10,2])
#apply(fdml$P, 3, eigen)

fdml$covar[,,1]

# perform cross validation 
dml.cv <- cvdmlpath(y,x,lambda = la)


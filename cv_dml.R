library(mvtnorm)

"%^%" =  function(x, n) 
  with(eigen(x), vectors %*% (values^n * t(vectors)))

# Fit model for a vector of lambda
dmlpath <- function(y,x,lambda){
  
  n <- dim(y)[1]      # nbr of observation 
  m <- dim(y)[2]     # nbr of response's component
  p <- dim(x)[2]     # nr of covariate 
  l.la <- length(lambda) # vector of lambda used for training
  
  # Initialize variable 
  coef <- array(0,dim=c(p,m,l.la))  
  covar <- array(0, dim = c(m,m,l.la))
  Phi <- array(0, dim = c(p,m,l.la))
  P <- array(0, dim = c(m,m,l.la))
  
  for (i in 1:l.la){
    fit <- JointSol(y,x,lambda=lambda[i])
    P[,,i] <- fit$P
    Phi[,,i] <- fit$Phi
    coef[,,i] <- fit$Phi %*% (fit$P %^% {-1})
    covar[,,i] <- fit$P %^% {-2}
  }
  
  #prepare results
  dimnames(coef) <- list(coefrow=1:p, coefcol = 1:m,lambda=lambda)
  dimnames(coef)[[2]] <- array('intercept',m)
  dimnames(Phi) <- list(Phirow = 1:p, Phicol = 1:m, lambda = lambda)
  dimnames(P) <- list(Prow = 1:m, Pcol =1:m, lambda = lambda)
  dimnames(covar) <- list(covrow = 1:m, covcol = 1:m, lambda= lambda) 

  res <- list(lambda=lambda,coef=coef,covar=covar, Phi =Phi, P = P)
  class(res) <- 'dmlpath'
  res
  
}

##Log Loss
logloss <- function(coef,covar,x,y){
  n <- dim(y)[1]
  dnregr <- vector()
  xcoef <- array(0,dim=c(n,m))
  for(i in 1:n){
      xcoef[i,] <- x[i,]%*%coef
      dnregr[i] <- dmvnorm(y[i,],mean = xcoef[i,], sigma = covar, log = FALSE) 
  }
  
  probdnregr <- dnregr
  dmix <- rowSums(probdnregr)
  -2*sum(log(dmix))
  
}


##Prediction Error 

predloss <- function(model,x,y){
  
  l <- length(model$lambda)
  coef <- model$coef
  covar <- model$covar
  loss <- rep(0,l)
  for (i in 1:l){
    loss[i] <- logloss((coef[,,i]),covar[,,i],x,y)
  }
  pred <- list(loss=loss,indexopt=which.min(loss))
  pred
}


### Cross Validation ###

cv.folds <- function (n, folds = 10) 
{
  split(sample(1:n), rep(1:folds, length = n))
}


cvdmlpath <- function(y,x,lambda,K=10){
  all.folds <- cv.folds(dim(y)[1],K)
  errmat <- matrix(0,length(lambda),K)
  for (i in seq(K)){
    omit <- all.folds[[i]]
    fit1 <- dmlpath(y=y[-omit,],x=x[-omit,],lambda=lambda)
    errmat[,i] <- predloss(fit1,x[omit,],y[omit,])$loss
    cat("\n CV Fold", i, "\n\n")
  }
  cv <- apply(errmat, 1, mean)
  cv.error <- sqrt(apply(errmat, 1, var)/K)
  dimnames(errmat) <- list(lambda=lambda,fold=1:K)
  
  object <- list(lambda = lambda, cv = cv, cv.error = cv.error,errmat=errmat)
  invisible(object)
}


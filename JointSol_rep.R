# Join update on the reparametrized value  

library(glasso)
library(lars)

JointSol = function(Y,X,lambda){
	
	N = dim(Y)[1]
	M = dim(Y)[2]
	P = dim(X)[2]

	## Centered Data

	cent.Y = matrix(NA,N,M)
	cent.X = matrix(NA,N,P)

	cent.Y = Y - matrix(rep(colMeans(Y),each=N),N,M)
	cent.X = X - matrix(rep(colMeans(X),each=N),N,P)
	
	## Estimating Phi and P
	update_PhiC = function(lambda){
	
	  ## Initial P
	  init.P = (t(cent.Y)%*%cent.Y)^{-0.5}/N
	  W = matrix(1,P,M)
	  
	  ## Initial Phi #TASK1

		init.Phi = matrix(0,P,M)
    init.YP = cent.Y %*% init.P
		for(i in 1:M){		
			lasso = lars(cent.X,init.YP[,i],use.Gram=FALSE,intercept=FALSE)
			init.Phi[,i] = predict.lars(lasso,cent.X,lambda/2,type="coefficients",mode="lambda")$coefficients
		}

		## Updating_PhiP
		output = JointLasso(cent.Y,cent.X,init.Phi,init.P,W,lambda)
		return(output)	
	}

	result = update_PhiC(lambda)
	resultPhi = result[["Phi"]]
	resultP = result[["P"]]
	
	return(list(Phi=resultPhi,P=resultP))
}

JointLasso = function(cent.Y,cent.X,init.Phi,init.P,W,lambda1){
  y<-cent.Y%*%init.P #1.should it be centralized?
  #2. no glasso for P initialization
  
	cur.Phi = init.Phi
	cur.P = init.P

	update = Onetime_update(cent.Y,cent.X,cur.Phi,cur.P,W,lambda1)
	updatePhi = update[["updatePhi"]]
	updateP = update[["updateP"]]
		
	maxDiff = max(max(abs(cur.Phi-updatePhi)),max(abs(cur.P-updateP)))
	iter_n = 1	

	while(maxDiff > 1.0e-3 && iter_n < 10)
	{ 
		cur.Phi = updatePhi
		cur.P = updateP
		
		update = Onetime_update(cent.Y,cent.X,cur.Phi,cur.P,W,lambda1)
		updatePhi = update[["updatePhi"]]
		updateP = update[["updateP"]]

		maxDiff = max(max(abs(cur.Phi-updatePhi)),max(abs(cur.P-updateP)))
		iter_n = iter_n+1
	}
	
	return(list(P = updateP, Phi = updatePhi, iteration = iter_n))	
}


Onetime_update = function(cent.Y,cent.X,cur.Phi,cur.P,W,lambda1){

	N = dim(cent.Y)[1]
	M = dim(cent.Y)[2]
	P = dim(cent.X)[2]
	dims = c(N,M,P)
  
	cent.YP<-cent.Y%*%cur.P #1.should it be centralized?
	
	## update P
  	A = t(cent.YP-cent.X%*%cur.Phi)%*%(cent.YP-cent.X%*%cur.Phi)/N
  	Cest = glasso(A,rho=1,penalize.diagonal=FALSE,wi.init=cur.P)
  	updateP = Cest$wi

	## Update Phi #TASK2
   
	updatePhi = numeric(P*M)
	iteration = numeric(1)
  
	dyn.load("UpdatePhi_rep.so")

	output = .C("UpdatePhi",as.double(t(cent.YP)),as.double(t(cent.X)),as.double(t(cur.Phi)),as.integer(dims),
	as.double(t(W)),as.double(lambda1),updatePhi=as.double(updatePhi),iteration=as.integer(iteration))

	dyn.unload("UpdatePhi_rep.so")

	updatePhi = matrix(output[["updatePhi"]],P,M,byrow=T)

	return(list(updateP=updateP,updatePhi=updatePhi))

}


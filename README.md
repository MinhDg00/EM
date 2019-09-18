
#### Research project supervised by Prof. Sunyoung Shin.
  Implement coordinate-descent algorithm for estimation and variable selection in reparametrized penalized Gaussian models. 

#### Main idea of the algorithm: 
Make effort to provide a efficient algorithm for numerical optimization with reparametrization and other convergent properties. The algorithm's purpose is to tackle with sparsity and large number of features  

#### File Instruction 
  * dmlsim.R  is the simulation script
  * JointSol_rep.R is revised DML algorithm by Lee & Liu to fit our case  
  * UpdatePhi_rep.c is the script to update the new parameter Phi 
  * cv_dml.R is the script containing cross validation as well as dmlpath function ( fit model for a vector of lambdas)
  * mvnorm.R is an open package for multivariate normal distribution 
  
#### Note
  * To run R code, you first need to compile the C code on your local machine 
  * With slightly large size of data, the code can run very slowly. The computation can be unstable for large number of features/ dimensions
  
#### Reference
```
Städler, Nicolas & Bühlmann, Peter & van de Geer, Sara. (2010). L1-Penalization for Mixture Regression Models. TEST: An Official Journal of the Spanish Society of Statistics and Operations Research. 19. 209-256. 10.1007/s11749-010-0197-z.
```
```
Lee, Wonyul & Liu, Yufeng. (2012). Simultaneous Multiple Response Regression and Inverse Covariance Matrix Estimation via Penalized Gaussian Maximum Likelihood. Journal of multivariate analysis. 111. 241-255. 10.1016/j.jmva.2012.03.013. 
```
####

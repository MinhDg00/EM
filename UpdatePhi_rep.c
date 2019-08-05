
# include <R.h>
# include <Rmath.h>
void UpdateOnephi(double *YP, double *X, double *Phi, int *dims, 
double *W, double *lambda, int *idx);

void UpdatePhi(double *YP, double *X, double *Phi, int *dims,
double *W, double *lambda, double *updatePhi, int *iteration)
{
    int i=1, j;
    int M = dims[1];
    int P = dims[2];
    int idx[2];
    int iter_n = 0;
    double maxDiff = 1;
    int iter_count = 0;
    
    while(maxDiff>1e-4 && iter_n<1000)
    {
        maxDiff = 0;
        
        for(i=0;i<M;i++)
        {
          for(j=0;j<P;j++)
          {
            idx[0] = i;
            idx[1] = j;

            double old_Phi = Phi[j*M+i];
            double tempDiff;

            UpdateOnephi(YP, X, Phi,dims,W,lambda,idx);
	    iter_count++;
            tempDiff = sign(old_Phi-Phi[j*M+i])*(old_Phi-Phi[j*M+i]);
            maxDiff = (maxDiff>tempDiff)? maxDiff:tempDiff;
          }
        }
      iter_n++;
    }
    
    for(i=0;i<M;i++)
      for(j=0;j<P;j++)
        updatePhi[j*M+i] = Phi[j*M+i];
    *iteration = iter_count;
}



void UpdateOnephi(double *YP, double *X, double *Phi, int *dims, 
double *W, double *lambda, int *idx)
{
    int i, j, k;
    int N = dims[0];
    int M = dims[1];
    int P = dims[2];
    int m = idx[0];
    int p = idx[1];
    double *residual = (double *)malloc(N*M*sizeof(double));
    double *presidual = (double *)malloc(M*sizeof(double));
    double sign_factor1,sign_factor2,sign_factor, updated_Phi;
    
    for(i=0;i<N;i++)
        for(j=0;j<M;j++)
        {
          residual[i*M+j]=0;
          for(k=0;k<P;k++)
              residual[i*M+j]+=X[i*P+k]*Phi[k*M+j];
	  residual[i*M+j]=YP[i*M+j]-residual[i*M+j];
        }
    for(i=0;i<M;i++)
    {
	     presidual[i]=0;
	     for(k=0;k<N;k++)
		      presidual[i]+=X[k*P+p]*residual[k*M+i];
    }

    for(i=0,sign_factor1=0;i<M;i++)
    	 sign_factor1+=presidual[i];

    for(i=0,sign_factor2=0;i<N;i++)
	     sign_factor2+=X[i*P+p]*X[i*P+p];
    
    sign_factor = sign_factor1/sign_factor2+Phi[p*M+m];

    updated_Phi = sign(sign_factor)*sign_factor-(*lambda)*W[p*M+m]/(2*sign_factor2);
    Phi[p*M+m] = sign(sign_factor)*updated_Phi*(updated_Phi>0);
    
    free(residual);
    free(presidual);
}

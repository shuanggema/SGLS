#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "mex.h"

int Sum(int *G, int n)
{
  int i,value=0;
  for( i = 0; i < n; i++)
     value += G[i];
     
  return value;
}

double Norm(double *z,int n)
{
  int i;
  double sum=0,value;
  for ( i = 0; i < n; i++)
  {
    sum += z[i]*z[i];
  }
  value = pow(sum,0.5);
  return value;
}

double Norm_b(double *Beta, int *G, int pos_s, int j)
{
  double sum = 0, value;
  int jj;
  for ( jj = 0; jj< G[j]; jj++ )
  {
     sum += Beta[pos_s+jj]*Beta[pos_s+jj];
  }
  value = pow(sum,0.5);
  return value;
}

double Inner_prod(double *X,double *Beta,int totalp)
{
  int j;
  double sum=0;
  
  for (j = 0; j < totalp; j++)
  {
    sum += X[j]*Beta[j];
  }
  return sum;
}

void S_func(double *z,double *X,double *Beta,double *r,double a, double c,
            int pos_s,int pos_e,int d_j,int n,int totalp)
{
  int i,j;
  double norm, cons, *diff,sum;//diff[d_j]
 
  diff = malloc(sizeof(double) * (int)d_j);
  norm = Norm(z,d_j);
  cons = (1 - c/norm);
  
  if (cons < 0)
    cons=0;
  
  for ( j = pos_s; j <= pos_e; j++)
  {
    diff[j-pos_s] = cons*z[j-pos_s]/a - Beta[j]; 
    Beta[j] = cons*z[j-pos_s]/a;
  }
  //update r here;
  for ( i = 0; i < n; i++)
  {
    sum = 0;
    for ( j = pos_s; j <= pos_e; j++)
    {
      sum += X[totalp*i+j]*diff[j-pos_s];
    }
    r[i] -= sum;
  }
}

void Update_beta(double *X, double *r, double *A, double *Beta, double *Lambda, 
                 int n,int totalp,int *pos_s,int *pos_e,int *G, int J, int P, int D)
{
  int i,j,k,s=pos_s[J],e=pos_e[J];
  double *z, a,c, sum1=0,sum2=0;//z[G[J]]
  
  z = malloc(sizeof(double) * G[J]);
  for ( j = s;j <= e;j++)
  {
    //let i-pos_s be the position for z[j];
    z[j-s] = 0;
    for ( i = 0; i < n; i++)
    {
      z[j-s] += X[totalp*i+j]*r[i]/n;
    }
    z[j-s] += Beta[j];
  }

  for ( k = 0; k < P; k++)
  {  
    if (  A[k*P+J]!=0 )
      sum1 += A[k*P+J];
  }  

  for ( k = 0 ; k < P; k++)
  {
    if (  A[k*P+J]!=0 )
      sum2 += A[k*P+J]*Norm_b(Beta,G,pos_s[k],k)/pow(G[k],0.5);
  }

  a = 1 + Lambda[1]*D/G[J]*sum1 - 1/Lambda[2];
  c = pow(G[J],0.5)*Lambda[0]-Lambda[1]*D/pow(G[J],0.5)*sum2;   
    
  S_func(z,X,Beta,r,a,c,s,e,G[J],n,totalp);
  
  if ( Norm_b(Beta,G,pos_s[J],J) > Lambda[0]*Lambda[2] ) //Lambda[0]:lambda1;Lambda[1]:lambda2;Lambda[2]:gamma;
  {
    a = 1 + Lambda[1]*D/G[J]*sum1;
    c = -Lambda[1]*D/pow(G[J],0.5)*sum2;
    S_func(z,X,Beta,r,a,c,s,e,G[J],n,totalp);
  } 
}

double gmcp_lapl(double *Y, double *X, int *G, double *A, int *Param, double *Lambda, double *Epsilon,
              double *Beta_ini,double *Beta)
{
  //int i,j,pos_s[Param[1]],pos_e[Param[1]],count=0,n=Param[0],p=Param[1],totalp=Param[2],iter=Param[3],d=Param[4];
  int i,j,*pos_s, *pos_e, count=0, n=Param[0],p=Param[1],totalp=Param[2],iter=Param[3],d=Param[4];
  //p is the number of group, totalp is the number of covariates;
  //double **Xt, y_t[Param[0]],r[Param[0]],old_beta[Param[2]], diff, diff_beta[Param[2]];
  double **Xt, *y_t, *r, *old_beta, diff, *diff_beta,tmp;
  //int *active_set;
  
  pos_s = malloc(sizeof(int) * (int)Param[1]);
  pos_e = malloc(sizeof(int) * (int)Param[1]);
  y_t = malloc(sizeof(double) * (int)Param[0]);
  r   = malloc(sizeof(double) * (int)Param[0]);
  old_beta = malloc(sizeof(double) * (int)Param[2]);
  diff_beta= malloc(sizeof(double) * (int)Param[2]);
  //active_set= malloc(sizeof(double) * (int)Param[2]);
  Xt = malloc(sizeof(double *) * (int)Param[0]); 
  for ( i = 0; i < Param[0]; i++)   
  { 
     Xt[i] = malloc(sizeof(double) * (int)Param[2]); 
  }
    
  for ( i = 0; i < Param[0]; i++ ){
    for ( j = 0; j < Param[2]; j++ ){
       Xt[i][j] = *(X+Param[2]*i+j); 
    }
  }   
  
  for( j = 0; j < p; j++){
    if ( j == 0)
       pos_s[j]=0;
    else
       pos_s[j]=Sum(G,j);
    
    pos_e[j]=pos_s[j]+G[j]-1;
  }
  // do initialization in outer loop;
  for( i = 0; i < n; i++){
    y_t[i] = Inner_prod(Xt[i],Beta_ini,totalp);
    r[i] = Y[i] - y_t[i];
  }
  
  for ( j = 0; j < totalp; j++){
      Beta[j] = Beta_ini[j];
      old_beta[j] = Beta[j];
  //    active_set[j]=1;
  }  
  
  do
  {
  
    for ( j = 0; j < totalp; j++)
      old_beta[j] = Beta[j];    
    
    for ( j = 0; j < p; j++){
    //   if ( active_set[j]!=0 ){
         Update_beta(X, r, A,Beta, Lambda, n, totalp, pos_s, pos_e, G,j,p,d);
    //     tmp = Norm_b(Beta,G,pos_s[j],j);
    //     if ( tmp < 1e-10 )
    //       active_set[j]=0;}
    }

    for ( j = 0; j < totalp; j++)
      diff_beta[j]=old_beta[j]-Beta[j];
      
    diff = Norm(diff_beta,totalp);        
    count=count+1;
  }
  while ( count <= iter && diff > Epsilon[0]); 
  printf("iter=%d\n",count); 
  
  for ( i = 0 ; i < n; i++)
  {  
    free(Xt[i]);
  }  
  free(Xt);
  return 0;
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )    
{ 
    int *G, *Param; 
    double *Y,*X,*A,*Lambda,*Epsilon,*Beta_ini,*Beta_out; 
    size_t m,n; 
    int beta_p;
    
    /* Check for proper number of arguments */
    
    if (nrhs != 8) { 
	    mexErrMsgIdAndTxt( "MATLAB:yprime:invalidNumInputs",
                "Two input arguments required."); 
    } 
    /*else if (nlhs > 1) {
	    mexErrMsgIdAndTxt( "MATLAB:yprime:maxlhs",
                "Too many output arguments."); 
    } */
    
    //Beta_OUT = mxCreateDoubleMatrix( (mwSize)m, (mwSize)n, mxREAL); 

    /* Assign pointers to the various parameters */ 
    Y = mxGetPr(prhs[0]); X = mxGetPr(prhs[1]);  G = (int * ) mxGetData(prhs[2]);
    A = mxGetPr(prhs[3]); Param = (int * ) mxGetData(prhs[4]); Lambda = mxGetPr(prhs[5]);
    Epsilon = mxGetPr(prhs[6]); Beta_ini = mxGetPr(prhs[7]);
    
    //beta_p=mxGetN(prhs[7]);
	beta_p=Param[2];
    plhs[0] = mxCreateDoubleMatrix( beta_p, 1, mxREAL); 
    Beta_out = mxGetPr(plhs[0]);
    /*Assign output pointers */
     
    /* Do the actual computations in a subroutine */
    //yprime(yp,t,y); 
    gmcp_lapl(Y, X, G, A, Param, Lambda, Epsilon, Beta_ini, Beta_out);
    //memcpy(Beta_out, Beta, sizeof(double)*beta_p);
    return;
    
}


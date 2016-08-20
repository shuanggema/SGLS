source("CallC.r")

lambda_max <- function(j,x,y,pos_s,pos_e,group)
{
  t=x[,pos_s[j]:pos_e[j]];
  n=nrow(x);
  inner=sum(colSums(t*as.vector(y)/n)^2)^0.5/group[j]^0.5
  list(inner=inner)
}

dyn.load("gmcp_lapl.so")

standard <- function (x) {
   p <- ncol(x)
   n <- nrow(x)
   x.mean <- matrix(rep(apply(x,2,mean),n),n,p,byrow=T)
   x.std <- t(t((x-x.mean))/(apply(x,2,sd)*(n-1)^0.5)*n^0.5)
}

gmcp_lapl <- function(x,y,beta,A,G,lambda1,lambda2,gamma,orthonormalize=F)
{
  n <- length(y)
  totalp <- dim(x)[2]
  p <- length(G)
  iter <- 200

  pos_e <- cumsum(G)
  pos_s <- pos_e - G + 1

  if ( orthonormalize == F)
  {
   fun <- function(i,x)
   {
     sigma = t(x[,pos_s[i]:pos_e[i]])%*%x[,pos_s[i]:pos_e[i]]/n
     R = chol(sigma)  
     x.t = x[,pos_s[i]:pos_e[i]]%*%solve(R)
     return(list(R=R,x.t=x.t))
   }
   result <- lapply(seq(1,p),fun,x) 

   fun <- function(j,result)
   {
     result[[j]]$x.t
   }
   x.t <- matrix(unlist(lapply(seq(1,p),fun,result)),nrow=n)
  }
  else if ( orthonormalize == T)
    x.t = x 

  beta <- numeric(totalp)
  param <- c(n, p, totalp,iter,max(G)) 
  epsilon <- 1E-6
  fit <- .C("GMCP_lapl", y=as.double(y),x=as.double(t(x.t)),G=as.integer(G),A=as.double(A),param=as.integer(param),
                         lambda=as.double(c(lambda1,lambda2,gamma)),epsilon=as.double(epsilon),beta=as.double(beta))
              
  beta <- numeric()
  b = fit$beta
  for ( j in 1:p )
  {
    if ( orthonormalize == F)
    {
      beta <- c(beta,solve(result[[j]]$R,b[pos_s[j]:pos_e[j]]))
    }
    else if (orthonormalize == T)
    {
      beta <- b
    }
  }
  list(beta=beta,b=b)
}

sp <- function(x,y,A,group,lambda2,gamma,ps,pe,n.step,epsilon)
{
  p <- length(group);
  n <- nrow(x);
  totalp <- ncol(x);

  y.c <- numeric(n)
  for ( i in 1:m)
    y.c[ps[i]:pe[i]] <- (y[ps[i]:pe[i]] - mean(y[ps[i]:pe[i]]));

  lambda.max=max(unlist(lapply(seq(1,p),lambda_max,x,y.c,pos_s,pos_e,group)));  

  lambda1.max <- lambda.max
  lambda1.min <- lambda1.max*epsilon
  ss <- (log(lambda1.max)-log(lambda1.min))/(n.step-1)
  lambda1.seq <- numeric(n.step)
  for ( i in 1:n.step )
  {
    lambda1.seq[i] <- exp(log(lambda1.max)-ss*(i-1))
  }

  beta.sp <- matrix(rep(0,n.step*totalp),nrow=n.step)
  for ( i in 1:n.step)
  {
    if ( i!=1 )
    {
     lambda1=lambda1.seq[i];
     fit=gmcp_lapl(x,y,beta.sp[i-1,],A,group,lambda1,lambda2,gamma,orthonormalize=T)
     beta.sp[i,]=fit$beta;
    }
    if ( i==1 )
    {
     lambda1=lambda1.seq[i];
     fit=gmcp_lapl(x,y,rep(0,totalp),A,group,lambda1,lambda2,gamma,orthonormalize=T)
     beta.sp[i,]=fit$beta; 
    }
    
  }
  list(beta.sp=beta.sp,lambda1.seq=lambda1.seq)
}


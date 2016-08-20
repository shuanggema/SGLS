

lambda_max <- function(j,x,y,pos_s,pos_e,group)
{
  t=x[,pos_s[j]:pos_e[j]];
  n=nrow(x);
  inner=sum(colSums(t*as.vector(y)/n)^2)^0.5/group[j]^0.5
  list(inner=inner)
}

fun <- function(i,x,y,A,m,fold.d,fold.g,n.fold,group,lambda1,lambda2,gamma)
{
  delete = unlist(fold.d[i])
  x.train = x[-delete,]
  y.train = y[-delete]
  x.valid = x[delete,]
  y.valid = y[delete]
  
  n.t <- numeric(m)
  for ( k in 1:m )
  {
    n.t[k] = 0
    for ( j in 1:n.fold)
    {
      if ( j != i )
        n.t[k] = n.t[k] + length(fold.g[[k]][[j]])
    }
  } 
  pe.t <- cumsum(n.t)
  ps.t <- pe.t - n.t + 1
  
  source("standard.r")
  index=seq(1,dim(x.train)[2])
  x.tnew=matrix(unlist(lapply(index,standardize,x.train)),nrow=dim(x.train)[1])
  
  y.tnew=numeric(length(y.train))
  for ( k in 1:m )
    y.tnew[ps.t[k]:pe.t[k]]=y.train[ps.t[k]:pe.t[k]]-mean(y.train[ps.t[k]:pe.t[k]])
  
  source("CallC.r")
  fit=gmcp_lapl(x.tnew,y.tnew,A,group,lambda1,lambda2,gamma,orthonormalize=T);  
  
  n.v <- numeric(m)
  for ( k in 1:m )
  {
    n.v[k] = 0
    for ( j in 1:n.fold)
    {
      if ( j == i )
        n.v[k] = n.v[k] + length(fold.g[[k]][[j]])
    }
  }
  pe.v <- cumsum(n.v)
  ps.v <- pe.v - n.v + 1
  
  y.c.valid <- numeric(length(y.valid))
  for ( i in 1:m)
    y.c.valid[ps.v[i]:pe.v[i]] <- y.valid[ps.v[i]:pe.v[i]] - mean(y.valid[ps.v[i]:pe.v[i]])
  
  index=seq(1,dim(x.valid)[2])
  x.vnew=matrix(unlist(lapply(index,standardize,x.valid)),nrow=dim(x.valid)[1])
  sse <- sum((y.c.valid - x.vnew%*%fit$beta)^2)
  list(sse=sse)
}

parCv <- function(cl,lambda1.seq,lambda2.seq,m,x,y,A,group,fold.d,fold.g,n.fold,gamma)
{
  step1 <- length(lambda1.seq)
  step2 <- length(lambda2.seq)
  PE <- matrix(rep(0,step1*step2),nrow=step1);
  index <- seq(1:n.fold)
  
  for( j in 1:step2)
  {
    lambda2=lambda2.seq[j]
    for ( i in 1: step1)
    {
      lambda1=lambda1.seq[i]
      v <- clusterApply(cl, index,fun,x,y,A,m,fold.d,fold.g,n.fold,group,lambda1,lambda2,gamma)
      PE[i,j] <- sum(unlist(v))/n.fold
    }
  }
  list(PE=PE)
}

cv.optim <- function(x,y,A,group,lambda2.seq,ps,pe,n.fold,epsilon,n.step,gamma,cl)
{
  p <- length(group);
  n <- nrow(x);
  totalp <- ncol(x);
  
  m <- length(ps);
  fold.g=vector("list",m);
  fold.d=vector("list",n.fold);
  for ( k in 1:m )
    fold.g[[k]] <- splitList(sample(ps[k]:pe[k]),n.fold);
  for ( k in 1:n.fold )
  {
    a <- numeric(0);
    for ( j in 1:m )
    {
      a <- c(a,fold.g[[j]][[k]]);
    }
    fold.d[[k]]=a;
  }
  
  y.c <- numeric(n)
  for ( i in 1:m)
    y.c[ps[i]:pe[i]] <- (y[ps[i]:pe[i]] - mean(y[ps[i]:pe[i]]));
  
  pos_e <- cumsum(group);
  pos_s <- pos_e - group + 1;
  d=rep(1,totalp);
  
  lambda.max=max(unlist(lapply(seq(1,p),lambda_max,x,y.c,pos_s,pos_e,group)));  
  
  lambda1.max <- lambda.max
  lambda1.min <- lambda1.max*epsilon
  ss <- (log(lambda1.max)-log(lambda1.min))/(n.step-1)
  lambda1.seq <- numeric(n.step)
  for ( i in 1:n.step )
  {
    lambda1.seq[i] <- exp(log(lambda1.max)-ss*(i-1))
  }
  
  
  cv = parCv(cl,lambda1.seq,lambda2.seq,m,x,y.c,A,group,fold.d,fold.g,n.fold,gamma)
  
  list(lambda1.seq=lambda1.seq,lambda2.seq=lambda2.seq,PE=cv$PE,fold.d=fold.d,fold.g=fold.g)
}


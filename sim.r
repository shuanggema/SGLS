match.fun.bal <- function(betaf,p)
{
  a=c(1:60)
  b=seq(1:p)
  b=b[betaf!=0]
  matched=match(a,b,nomatch=-1)
  tp <- length(matched[matched!=-1])
  fn <- length(a) -tp
  fp <- length(b) - tp
  tn <- p - tp -fn -fp
  list(tp=tp,fn=fn,fp=fp,tn=tn)
}

weighted.standard.x <- function (x,w) {
   p <- ncol(x)
   n <- nrow(x)
   x.wmean <- matrix(rep(apply(x,2,weighted.mean,w),n),n,p,byrow=T)
   x.std <- (x-x.wmean)*w^0.5
}
weighted.standard.y <- function (y,w) {
   n <- length(y)
   y.wmean <- weighted.mean(y,w)
   y.std <- (y-y.wmean)*w^0.5
}

transform <- function(x,y,w)
{
  x.sd=weighted.standard.x(x,w); 
  y.sd=weighted.standard.y(y,w);  
  list(x.sd=x.sd,y.sd=y.sd)
}

standard <- function (x) {
   p <- ncol(x)
   n <- nrow(x)
   x.mean <- matrix(rep(apply(x,2,mean),n),n,p,byrow=T)
   x.std <- t(t((x-x.mean))/(apply(x,2,sd)*(n-1)^0.5)*n^0.5)
}

standardize  <- function(x) (x-mean(x))/sd(x)
process <- function(xt,m,n)
{
  ntotal=nrow(xt)
  p=ncol(xt)
  pe=cumsum(n)
  ps=pe-n+1

  x=matrix(rep(0,ntotal*p*m),nrow=ntotal)
  for( j in 1:p )
    for ( k in 1:m )
    {
      start=ps[k]
      end = pe[k]
      x[start:end,m*(j-1)+k]=xt[start:end,j]
    }

  x0=matrix(rep(0,ntotal*m),nrow=ntotal)
  for ( k in 1:m )
  {
     start=ps[k]
     end = pe[k]
     x0[start:end,k]=rep(1,end-start+1)
  }
  list(x0=x0,x=x,ps=ps,pe=pe)
}

can.cor <- function(j,k,x,pos_s,pos_e)
{
  a1 <- x[,(pos_s[j]:pos_e[j])]
  a2 <- x[,(pos_s[k]:pos_e[k])]
  cc <- cancor(a1,a2)
  max(cc$cor)
}

vec.cor <- function(j,k,x,pos_s,pos_e,m)
{
  a1 <- x[,(pos_s[j]:pos_e[j])]
  a2 <- x[,(pos_s[k]:pos_e[k])]
  abs(cor(rowSums(a1)/m^0.5,rowSums(a2)/m^0.5))
}

N.1 <- function(x,p,m,pos_s,pos_e,c)
{
  A = matrix(rep(0,p*p),nrow=p)
  cutoff=(exp(2*c/(nrow(x)-3)^0.5)-1)/(exp(2*c/(nrow(x)-3)^0.5)+1)
  for ( j in 1:(p-1) )
    for ( k in (j+1):p )
    {
      tmp=vec.cor(j,k,x,pos_s,pos_e,m);
      if ( abs(tmp)>=cutoff )
        A[j,k]=1
    }

  for ( k in 1:(p-1) )
    for ( j in (k+1):p )
     A[j,k]=A[k,j];
  list( A=A)
}

N.6 <- function(x,p,m,pos_s,pos_e,c)
{
  A = matrix(rep(0,p*p),nrow=p)
  cutoff=(exp(2*c/(nrow(x)-3)^0.5)-1)/(exp(2*c/(nrow(x)-3)^0.5)+1)
  for ( j in 1:(p-1) )
    for ( k in (j+1):p )
    {
      tmp=vec.cor(j,k,x,pos_s,pos_e,m);
      if ( abs(tmp)>=cutoff )
        A[j,k]=abs(tmp)
    }

  for ( k in 1:(p-1) )
    for ( j in (k+1):p )
     A[j,k]=A[k,j];
  list( A=A)
}

N.2 <- function(x,p,m,pos_s,pos_e,cutoff)
{
  A = matrix(rep(0,p*p),nrow=p)
  
  for ( j in 1:(p-1) )
    for ( k in (j+1):p )
    {
      tmp=can.cor(j,k,x,pos_s,pos_e);
      if ( tmp>=cutoff )
        A[j,k]=1
    }

  for ( k in 1:(p-1) )
    for ( j in (k+1):p )
     A[j,k]=A[k,j];
  list( A=A)
}

N.7 <- function(x,p,m,pos_s,pos_e,cutoff)
{
  A = matrix(rep(0,p*p),nrow=p)
  
  for ( j in 1:(p-1) )
    for ( k in (j+1):p )
    {
      tmp=can.cor(j,k,x,pos_s,pos_e);
      if ( tmp>=cutoff )
        A[j,k]=tmp
    }

  for ( k in 1:(p-1) )
    for ( j in (k+1):p )
     A[j,k]=A[k,j];
  list( A=A)
}

N.3 <- function(x,p,m,pos_s,pos_e,alpha,tau0)
{
  A = matrix(rep(0,p*p),nrow=p)
  
  for ( j in 1:(p-1) )
    for ( k in (j+1):p )
    {
      tmp=can.cor(j,k,x,pos_s,pos_e);
      A[j,k]=1/(1+exp(-alpha*(tmp-tau0)));
    }

  for ( k in 1:(p-1) )
    for ( j in (k+1):p )
     A[j,k]=A[k,j];
  list( A=A)
}

N.4 <- function(x,p,m,pos_s,pos_e,alpha)
{
  A = matrix(rep(0,p*p),nrow=p)
  
  for ( j in 1:(p-1) )
    for ( k in (j+1):p )
    {
      tmp=can.cor(j,k,x,pos_s,pos_e);
      A[j,k]=tmp^alpha;
    }

  for ( k in 1:(p-1) )
    for ( j in (k+1):p )
     A[j,k]=A[k,j];
  list( A=A)
}

N.5 <- function(x,p,m,pos_s,pos_e,alpha,cutoff)
{
  A = matrix(rep(0,p*p),nrow=p)
  
  for ( j in 1:(p-1) )
    for ( k in (j+1):p )
    {
      tmp=can.cor(j,k,x,pos_s,pos_e);
      if ( tmp <= cutoff)
        tmp=0;
      A[j,k]=tmp^alpha;
    }

  for ( k in 1:(p-1) )
    for ( j in (k+1):p )
     A[j,k]=A[k,j];
  list( A=A)
}

locate.min <- function(X)
{
  opt.X=min(X);n=nrow(X);p=ncol(X);
  index=seq(1,n*p);
  indx=index[as.numeric(X)==opt.X];
  mod.X=indx%%n;
  if ( mod.X != 0 )
  {
     j=as.integer(indx/n)+1;
     i=mod.X;
  }
  else if ( mod.X == 0 )
  {
    j=as.integer(indx/n);
    i=n;
  }
  list(i=i,j=j)
}

data.process <- function(sdm,mag,m,n,p,group)
{
    pos_e = cumsum(group)
    pos_s = pos_e -group + 1

    data=generate(sdm,mag,n,p,0.9);
    x1=data$x1;y1=data$y1;d1=data$d1;w1=data$w1;
    x2=data$x2;y2=data$y2;d2=data$d2;w2=data$w2;
    x3=data$x3;y3=data$y3;d3=data$d3;w3=data$w3;

    w1.s <- w1[w1!=0];x1.s <- x1[w1!=0,];y1.s=y1[w1!=0];n1=length(w1.s);
    w2.s <- w2[w2!=0];x2.s <- x2[w2!=0,];y2.s=y2[w2!=0];n2=length(w2.s);
    w3.s <- w3[w3!=0];x3.s <- x3[w3!=0,];y3.s=y3[w3!=0];n3=length(w3.s);

    data1=transform(x1.s,log(y1.s),w1.s);
    data2=transform(x2.s,log(y2.s),w2.s);
    data3=transform(x3.s,log(y3.s),w3.s);
    
    x1.std=standard(data1$x.sd);
    x2.std=standard(data2$x.sd);
    x3.std=standard(data3$x.sd);

    y1.new=data1$y.sd-mean(data1$y.sd);
    y2.new=data2$y.sd-mean(data2$y.sd);
    y3.new=data3$y.sd-mean(data3$y.sd);
 
    xt=rbind(x1.std,x2.std,x3.std)
    y=c(y1.new,y2.new,y3.new)
    n.t=c(n1,n2,n3);
    data=process(xt,m,n.t)
    x=data$x; x0=data$x0; ps=data$ps; pe=data$pe
    totalp = ncol(x);

    source("standard.r")
    index=seq(1,totalp)
    x.new=matrix(unlist(lapply(index,standardize,x)),nrow=nrow(x))

    y.tmp=numeric(length(y))
    for ( k in 1:m )
      y.tmp[ps[k]:pe[k]]=y[ps[k]:pe[k]]-mean(y[ps[k]:pe[k]])

    y.new=y.tmp

 
   list(x.new=x.new,y.new=y.new,ps=ps,pe=pe)
}


sim <- function(n,p,sdm,mag,cl,n.sim,epsilon,gamma,N.type,example)
{
  m=3
  tp <- numeric()
  tn <- numeric()
  fp <- numeric()
  fn <- numeric()
  lambda1.opt.seq <- numeric()
  lambda2.opt.seq <- numeric()
  lambda.opt.seq <- numeric()
  lambda.min.seq <- numeric()
  lambda.max.seq <- numeric()

  group=rep(m,p)
  n.fold=5
  n.step=50

  for ( i in 1:n.sim )
  {
    data=data.process(sdm,mag,m,n,p,group);
    x.new=data$x.new;y.new=data$y.new; ps=data$ps;pe=data$pe; 

    pos_e = cumsum(group)
    pos_s = pos_e -group + 1
    
    if ( N.type == 1)
      A=N.1(x.new,p,m,pos_s,pos_e,3.29)$A;
    if ( N.type == 2)
      A=N.2(x.new,p,m,pos_s,pos_e,0.43002431)$A;  
    if ( N.type == 3)
      A=N.3(x.new,p,m,pos_s,pos_e,10,0.43002431)$A;
    if ( N.type == 4)
      A=N.4(x.new,p,m,pos_s,pos_e,6)$A;
    if ( N.type == 5)
      A=N.5(x.new,p,m,pos_s,pos_e,6,0.43002431)$A;
    if ( N.type == 6)
      A=N.6(x.new,p,m,pos_s,pos_e,3.29)$A;
    if ( N.type == 7)
      A=N.7(x.new,p,m,pos_s,pos_e,0.43002431)$A;

    source("cv.r")
    lambda2.seq=c(0,0.001,0.01,0.1,1,10);
    
    repeat
    {
     cv =cv.optim(x.new,y.new,A,group,lambda2.seq,ps,pe,n.fold,epsilon,n.step,gamma,cl)
     lambda1.seq=cv$lambda1.seq;
     lambda2.seq=cv$lambda2.seq;
     loc=locate.min(cv$PE);
    
     lambda1.opt.seq <- c(lambda1.opt.seq,lambda1.seq[loc$i])
     lambda2.opt.seq <- c(lambda2.opt.seq,lambda2.seq[loc$j])
    
     source("CallC.r")
     fit=gmcp_lapl(x.new,y.new,A,group,lambda1.seq[loc$i],lambda2.seq[loc$j],gamma,orthonormalize=T);
     PE.g=cv$PE[,1];  
     fit.g=gmcp_lapl(x.new,y.new,A,group,lambda1.seq[PE.g==min(PE.g)],0,gamma,orthonormalize=T);

     a=match.fun.bal(fit$beta,m*p)
     a.g=match.fun.bal(fit.g$beta,m*p)
     if ( a$fp/3<=40)
       break
    }
    #prediction on another sample with the same size and censoring scheme;
    data1=data.process(sdm,mag,m,n,p,group);
    x.new=data1$x.new;y.new=data1$y.new;
    mspe1=sum((y.new-x.new%*%fit$beta)^2);
    mspe2=sum((y.new-x.new%*%fit.g$beta)^2);
    
    lambda.opt.seq <- c(lambda.opt.seq,lambda1.seq[PE.g==min(PE.g)])
    lambda.min.seq <- c(lambda.min.seq,min(lambda1.seq))
    lambda.max.seq <- c(lambda.max.seq,max(lambda1.seq))
    
    
    index=seq(1,m*p)
    #group.s = group.id(beta,group)$group.s
    outfile1 <- gsub("( )", "", paste("gene_s_GMCP_lapl_",n,"_",p,"_",gamma,"_",example,".txt"))
    outfile2<- gsub("( )", "", paste("gene_s_GMCP_",n,"_",p,"_",gamma,"_",example,".txt"))
    outfile3 <- gsub("( )", "", paste("MSPE_GMCP_lapl_",n,"_",p,"_",gamma,"_",example,".txt"))
    outfile4<- gsub("( )", "", paste("MSPE_GMCP_",n,"_",p,"_",gamma,"_",example,".txt"))
    write.table(matrix(index[fit$beta!=0],nrow=1),outfile1,sep=" ",quote=F, row.names=F,col.names=F,append=T)
    write.table(matrix(index[fit.g$beta!=0],nrow=1),outfile2,sep=" ",quote=F, row.names=F,col.names=F,append=T)
    write.table(mspe1,outfile3,sep=" ",quote=F, row.names=F,col.names=F,append=T)
    write.table(mspe2,outfile4,sep=" ",quote=F, row.names=F,col.names=F,append=T)
    
    tp = c(tp,a$tp,a.g$tp);tn = c(tn,a$tn,a.g$tn);fp = c(fp,a$fp,a.g$fp);fn = c(fn,a$fn,a.g$fn);
  }
  list(tp=tp,tn=tn,fp=fp,fn=fn,lambda1.opt.seq=lambda1.opt.seq,lambda2.opt.seq=lambda2.opt.seq,
      lambda.opt.seq=lambda.opt.seq,lambda.min.seq=lambda.min.seq,lambda.max.seq=lambda.max.seq)
}



library(mvtnorm)
generate <- function(k,mag,n,p,quan)
{
 beta1 <- c(rep(0.5,5),rep(0,p-5))
 beta2 <- c(rep(0.5,5),rep(0,p-5))
 beta3 <- c(rep(0.5,5),rep(0,p-5))
 beta0 <- 0.5
 
 clusters <- as.integer(p/5)
 covar <- matrix(rep(0,p*p),nrow=p)
 
 for ( clus in 1:clusters)
 {
   for ( i in (5*clus-4):(5*clus))
     for ( j in (5*clus-4):(5*clus))
     covar[i,j]=mag^(abs(i-j))
 }

 x1 <- rmvnorm(n = n, mean = rep(0, p), covar) 
 x2 <- rmvnorm(n = n, mean = rep(0, p), covar)
 x3 <- rmvnorm(n = n, mean = rep(0, p), covar)
 
 T1=beta0 + x1%*%beta1 + rnorm(n = n, 0 , k)
 T2=beta0 + x2%*%beta2 + rnorm(n = n, 0 , k)
 T3=beta0 + x3%*%beta3 + rnorm(n = n, 0 , k)
 
 C1 <- runif(n,min=0,max=quantile(exp(T1),c(quan)))
 y1 <- pmin(exp(T1), C1)
 d1 <- as.numeric(exp(T1) <= C1)
 xs1=x1[order(y1),]
 d1=d1[order(y1)]
 y1=sort(y1)

 C2 <- runif(n,min=0,max=quantile(exp(T2),c(quan)))
 y2 <- pmin(exp(T2), C2)
 d2 <- as.numeric(exp(T2) <= C2)
 xs2=x2[order(y2),]
 d2=d2[order(y2)]
 y2=sort(y2)
 
 C3 <- runif(n,min=0,max=quantile(exp(T3),c(quan)))
 y3 <- pmin(exp(T3), C3)
 d3 <- as.numeric(exp(T3) <= C3)
 xs3=x3[order(y3),]
 d3=d3[order(y3)]
 y3=sort(y3)
 
 w1 <- numeric(n)
 w1[1]=d1[1]/n
 for ( i in 2:n )
 {
   tmp = 1
   for ( j in 1: (i-1) )
     tmp = tmp*((n-j)/(n-j+1))^d1[j]
     
   w1[i]=d1[i]/(n-i+1)*tmp
 }

 w2 <- numeric(n)
 w2[1]=d2[1]/n
 for ( i in 2:n )
 {
   tmp = 1
   for ( j in 1: (i-1) )
     tmp = tmp*((n-j)/(n-j+1))^d2[j]

   w2[i]=d2[i]/(n-i+1)*tmp
 }

 w3 <- numeric(n)
 w3[1]=d3[1]/n
 for ( i in 2:n )
 {
   tmp = 1
   for ( j in 1: (i-1) )
     tmp = tmp*((n-j)/(n-j+1))^d3[j]

   w3[i]=d3[i]/(n-i+1)*tmp
 }

 list(y1=y1,y2=y2,y3=y3,x1=xs1,x2=xs2,x3=xs3,
      d1=d1,d2=d2,d3=d3,w1=w1,w2=w2,w3=w3)
}


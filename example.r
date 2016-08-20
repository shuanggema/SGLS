source("generate_sp.r")

n=100
p=30
m=3
group=rep(m,p)
pos_e = cumsum(group)
pos_s = pos_e -group + 1

source("sim.r")

set.seed(10001)
data=data.process(2,0.9,m,n,p,group);
x.new=data$x.new;y.new=data$y.new; ps=data$ps;pe=data$pe; 
A=N.5(x.new,p,m,pos_s,pos_e,6,0.43002431)$A;

source("sp.r")

sp.sgls=sp(x.new,y.new,A,group,lambda2=0.1,100000,ps,pe,50,0.2)
b=sp.sgls$beta.sp

source("CallC.r")
fit.l=gmcp_lapl(x.new,y.new,A,group,0.01712265,0,6,orthonormalize=T)

######################################################################

png("sgls.png",width=480,height=480)
plot(sp.sgls$lambda1.seq,b[,1],ylim=c(min(b),max(b)),lty=1,type = 'l',col="blue",axes=F,ann=FALSE)
for ( i in c(2:3) )
 points(sp.sgls$lambda1.seq,b[,i],col="blue",lty=1,type = 'l')
for ( i in c(4:6) )
 points(sp.sgls$lambda1.seq,b[,i],col="red",lty=1,type = 'l')
for ( i in c(7:9) )
 points(sp.sgls$lambda1.seq,b[,i],col="yellow",lty=1,type = 'l')
for ( i in c(10:12) )
 points(sp.sgls$lambda1.seq,b[,i],col="orange",lty=1,type = 'l')
for ( i in c(13:15) )
 points(sp.sgls$lambda1.seq,b[,i],col="black",lty=1,type = 'l')
for ( i in c(16:90) )
 points(sp.sgls$lambda1.seq,b[,i],lty=3,type = 'l')
axis(2, at=seq(-0.01,0.03,0.01), line=0,las=2,cex.axis=1.5)#
axis(1, at=seq(0.02,0.08,0.02),line=0,las=1,cex.axis=1.5)#
mtext(expression(paste(beta)),side=2,line=3,las=2,cex=2)
title(xlab=expression(paste(lambda)),cex.lab =2)
dev.off()


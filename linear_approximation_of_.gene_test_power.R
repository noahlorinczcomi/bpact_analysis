m=100
ar1=function(n,rho) rho^toeplitz(0:(n-1))
tr=function(x) sum(diag(x))
R=ar1(m,0.5)
e0=m
v0=2*tr(R%*%R)
alpha0=e0^2/v0
beta0=e0/v0
alpha=0.05
Q=qgamma(1-alpha,shape=alpha0,rate=beta0)
mus=seq(0,2,length.out=100)
e1s=e0+m*mus^2
v1s=c()
for(i in 1:length(mus)) {
  mui=rep(mus[i],m)
  v1s[i]=v0+4*t(mui)%*%R%*%mui
}
alpha1=e1s^2/v1s
beta1=e1s/v1s
# plotting example of power
curve(dgamma(x,shape=alpha0,rate=beta0),50,300,xlab=expression('T'[k]),ylab='density')
s=seq(Q,qgamma(1-1e-10,shape=alpha0,rate=beta0),length.out=30)
y=dgamma(s,shape=alpha0,rate=beta0)
polygon(c(s,rev(s)),c(rep(0,length(s)),rev(y)),col='#708DFF5C',border='#708DFF5C')
curve(dgamma(x,shape=alpha0,rate=beta0),50,300,add=T)
k=50 # index to choose for example
curve(dgamma(x,shape=alpha1[k],rate=beta1[k]),50,400,add=T,lty=3)
s=seq(Q,qgamma(1-1e-10,shape=alpha1[k],rate=beta1[k]),length.out=30)
y=dgamma(s,shape=alpha1[k],rate=beta1[k])
polygon(c(s,rev(s)),c(rep(0,length(s)),rev(y)),col='#FF70705C',border='#FF70705C')
# plotting piecewise approximation to power
gamma=pgamma(Q,shape=alpha1,rate=beta1,lower.tail=FALSE)
plot(mus,gamma,type='l',lwd=2,xlab='ncp',ylab='power')
lines(c(0,0.3),c(alpha,0.1),col='red'); points(0.3,0.1,col='red',pch=19)
lines(c(0.3,0.9),c(0.1,0.975),col='red'); points(0.9,0.975,col='red',pch=19)
lines(c(0.9,2),c(0.975,1),col='red'); points(2,1,col='red',pch=19)

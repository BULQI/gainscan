library(minpack.lm)

g<-1.5
a<-1
d<-100
b<-2
x<-c(1:5000)
c<-600
y<-d+(a-d)/(1+(x/c)^b)^g

plot(c(1,5000), c(1,100), log="x")

lines(x,y, col=1)

##
c<-200
set.seed(20)
y<-d+(a-d)/(1+(x/c)^b)^g+rnorm(length(x),0, 0.01)

lines(x,y, col=3)

nls(y~d+(a-d)/(1+(x/c)^b)^g, start=list(a=1.1, b=2.1, c=200.5,d=90, g=1.0), trace=T)
gnls(y~d+(a-d)/(1+(x/c)^b)^g, start=list(a=9.1, b=4.1, c=200.5,d=10, g=1.0),verbose=T)

####mle estimation

LL<-function(a0,b0,c0,d0,g0, sigma0)
{
	r<- y-(d0+(a0-d0)/(1+(x/c0)^b0)^g0)
	r<-dnorm(r,0, sigma0, log=TRUE)
	-sum(r)
}

LL(a0=1.1,b0=2.1,c0=200.5, d0=99, g0=1.5,sigma0=0.1)

fit<-mle2(LL, start=list(a0=1.0, b0=2.0, c0=200,d0=100, g0=1.5, sigma0=0.01),fixed=list(a0=1,d0=100,b0=2,g0=1.5))#, nobs=length(y))


fit<-nls(y~SSlogis(x,Asym, xmid,scal))


#################now fit a different fpl
#y<-a+(d-a)/(1+exp((xmid-input)/scal))
#y<-a+(d-a)/(1+exp((xmid-input)/scal))^g

g<-2
a<-1
d<-100
xmid<-500
scal<- 50
xa<-seq(0,1000,by=10)

ya<-a+(d-a)/(1+exp((xmid-xa)/scal))^g#+rnorm(length(x),0, 0.1)

fitMin<-nlsLM(ya~a+(d-a)/(1+exp((xmid-xa)/scal)), start=list(a=2,d=20, xmid=100, scal=29))
###-------derivative of four pL
#dydx<-(d-a)/scal*(1+exp((xmid-xa)/scal))^(-2)*exp((xmid-xa)/scal)
dydx<-(d-a)/scal*(1/(1+exp((xmid-xa)/scal )))*(1-(1/(1+exp((xmid-xa)/scal ))))
dy2dx<-dydx/scal*(1-2*(1/(1+exp((xmid-xa)/scal ))))

###derivative of five pl
dydx5pl<-(d-a)*(g/scal)*(1/(1+exp((xmid-xa)/scal)))^g *(1-(1/(1+exp((xmid-xa)/scal))))
dy2dx5pl<-dydx5pl/scal*(g-(g+1)*(1/(1+exp((xmid-xa)/scal ))))

plot(c(-500,2000),c(-4e-3,4e-3),col=1, log="", type="n")
lines(xa, ya, col=2)
lines(xa,dydx,col=3, lty=2)
lines(c(xmid,xmid),c(0.01,100), lty=3, col=3)
lines(xa,dy2dx,col=4, lty=1)
lines(xa,dydx5pl, col=5,lty=2)
lines(xa,dy2dx5pl,col=6,lty=3)

######4pl with different model
#y<-a+(d-a)/(1+(x/c)^b)



g<-2
a<-1
d<-100
c<-500
b<- 50
xc<-seq(0,2000,by=50)

yc<-d+(a-d)/(1+(xc/c)^b)

dycdxc<--1*(a-d)*(b/c)*(1/(1+(xc/c)^b))*(1-(1/(1+(xc/c)^b)))*(1/(xc/c))

plot(xc,yc)


fit2<-nls(ya~SSfpl(xa,a,d,xmid, scal))


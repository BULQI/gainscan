library(ARPPA)

nGene<-1000
nTreatment<-2   #number of different beta
sampleSize<-20 ##this is the number of repeats for each group

alpha.mean<-0  #variance for alpha prior
alpha.sigma<-5
beta.mean<-0
beta.sigma<-5  #variance for beta prior

#v0<-10  # the unscaled factor for gamma given gamma<=>0
gammaN0.sigma<-10

p_diff<-0.01

#priors for variance distribution
d0<-4
s0<-2

#call it
set.seed(2004);
dataExp<-simulateExpression(nGene, nTreatment, sampleSize,
					alpha.mean=alpha.mean, alpha.sigma=alpha.sigma,
					beta.mean=beta.mean, beta.sigma=beta.sigma,
					prob.nonZero=0.01, gamma.sigma=gammaN0.sigma,
					epsilon.d0=d0, epsilon.s0=s0
					)

##now start fitting the model to estimate the EBayes parameters
#defining a S3 function for calculating the prior variance and df 
#based on the observed variance
calculatePriorVariance<-function (x,df)
{
	eg<-log(x)-digamma(df/2)+log(df/2)
	eg_bar<-mean(eg)
	tg_d02<-mean((eg-eg_bar)*(eg-eg_bar)*length(x)/(length(x)-1)-trigamma(df/2))
	d0<-trigammaInverse(tg_d02)*2
	s02<-exp(eg_bar+digamma(d0/2)-log(d0/2))
	
	prior<-list(d0=d0, s02=s02)
	
} 

scaledChiSq_Pdf<-function(x, d0, s0)
{
	x_t<-d0*s0*s0/(x)
	#x_t<-sqrt(x_t)
	1/(2^(d0/2)*gamma(d0/2))*x_t^(d0/2-1)*exp(-1*x_t/2)
	
	#dchisq(x_t, d0)
}
scaledChiSq_Pdfr<-function(x, d0, s0)
{
	x_t<-d0*s0*s0/(x)
	#x_t<-sqrt(x_t)
	#1/(2^(d0/2)*gamma(d0/2))*x_t^(d0/2-1)*exp(-1*x_t/2)
	
	dchisq(x_t, d0)
}
scaledChiSq_PdfGamma<-function(x, d0, s0)
{
	#x_t<-d0*s0*s0/(x)
	#x_t<-sqrt(x_t)
	#1/((2/(d0*s0*s0))^(d0/2)*gamma(d0/2))*x^(-1*(d0/2-1))*exp(-1*(x^-1/(2/(d0*s0*s0))))
	
	1/((2/(d0*s0*s0))^(d0/2)*gamma(d0/2))*(x)^((d0/2-1))*exp(-1*(x/(2/(d0*s0*s0))))
	
	#dchisq(x_t, d0)
}
scaledChiSq_PdfGammaR<-function(x, d0, s0)
{
	#x_t<-d0*s0*s0/(x)
	#x_t<-sqrt(x_t)
	#1/((2/(d0*s0*s0))^(d0/2)*gamma(d0/2))*x^(-1*(d0/2-1))*exp(-1*(x^-1/(2/(d0*s0*s0))))
	
	#1/((2/(d0*s0*s0))^(d0/2)*gamma(d0/2))*(x)^((d0/2-1))*exp(-1*(x/(2/(d0*s0*s0))))
	dgamma(x,shape=d0/2 ,scale=2/(d0*s0*s0))
	#dchisq(x_t, d0)
}

rGamma<-function(n, d0, s0)
{
	#x_t<-d0*s0*s0/(x)
	#x_t<-sqrt(x_t)
	#1/((2/(d0*s0*s0))^(d0/2)*gamma(d0/2))*x^(-1*(d0/2-1))*exp(-1*(x^-1/(2/(d0*s0*s0))))
	
	#1/((2/(d0*s0*s0))^(d0/2)*gamma(d0/2))*(x)^((d0/2-1))*exp(-1*(x/(2/(d0*s0*s0))))
	rgamma(n,shape=d0/2 ,scale=2/(d0*s0*s0))
	#dchisq(x_t, d0)
}
invGamma_Pdf<-function(x, d0, s0)
{
	#x_t<-d0*s0*s0/(x)
	#x_t<-sqrt(x_t)
	#1/((2/(d0*s0*s0))^(d0/2)*gamma(d0/2))*x^(-1*(d0/2-1))*exp(-1*(x^-1/(2/(d0*s0*s0))))
	
	#1/((2/(d0*s0*s0))^(d0/2)*gamma(d0/2))*(x)^((d0/2-1))*exp(-1*(x/(2/(d0*s0*s0))))
	dinvgamma(x,alpha=d0/2 ,beta=(d0*s0*s0)/2)
	#dchisq(x_t, d0)
}

####testing code
###following the example for R help for eBayes
#  See also lmFit examples

Y_ijk<-dataExp$exp
tY<-data.frame(t(Y_ijk))
tY$group<-c(rep(1,sampleSize),rep(2,sampleSize))
varY<-aggregate(tY,by=list(tY$group),FUN=var) #by doing aggregate, the new dataframe is arranged to have first column for "by" criteria.
		#that is why in the next step, we need to get rid of the first column and last one, since the last one is the one we added as the group info.
varK<-varY[,c(-1,-nGene-2)]

x<-as.numeric(c(varK[1,],varK[2,]))
df<-rep(sampleSize-1,length(x))

r2<-calculatePriorVariance(x,df)

xx<-seq(0,60, by=0.1)

y<-scaledChiSq_Pdf(xx,d0, s0 )
yr<-scaledChiSq_Pdfr(xx,d0, s0 )
#plot(density((x)))
#lines((xx), y, col=2, lty=2)

yy<-rScaledChisq(1000000, d0, s0*s0)
yyg<-1/rGamma(1000000, d0, s0)
yg<-scaledChiSq_PdfGamma(1/xx,d0,s0)
ygr<-scaledChiSq_PdfGammaR(1/xx,d0,s0)
yig<-invGamma_Pdf(xx, d0, s0)
plot(c(0,60), c(0, 0.1),type="n", xlab="var", ylab="density")
hist(yy,breaks=100,prob=TRUE, xlim=c(0,60))
lines(density(yy))
lines(density(yyg), col=2, lwd=2)
lines((xx), y, col=2, lty=3)
lines((xx),yr,col=3, lty=2)
lines((xx), yig, col=4, lty=3, lwd=2)
#lines(xx,yg, col=3, lty=4, lwd=2)
lines(xx,ygr, col=2, lty=4, lwd=2)
#plot((xx), y, col=2, lty=2)

chisq_pdf<-function(x, d0)
{
	1/(2^(d0/2)*gamma(d0/2))*x^(d0/2-1)*exp(-1*x/2)

}
ychi<-rchisq(10000,2)
plot(density(ychi))
xxchi<-seq(1,50,by=0.1)
yychi<-chisq_pdf(xxchi, 2)
lines(xxchi,yychi, lty=2, col=2)


######################code to run normal distribution fitting density

 xseq<-seq(-10,10,.01)
 y<- rnorm(10000,2,5.5)
  
 hist(y, prob=TRUE,  breaks=20)#ylim=c(0,.06),)
 curve(dnorm(x, mean(y), sd(y)), add=TRUE, col="darkblue", lwd=2)
 curve(dnorm(x, 2, 5.5), add=TRUE, col="darkred", lwd=2)

 lines(xseq, dnorm(xseq, 2,5.5), col=2,lty=2)
 lines(density(y), col=3)
 lines(xseq, dnorm(xseq, 0, sqrt(5.5)), col=2)



normal_Pdf<-function(x, m,s2)
{
	1/sqrt(2*pi*s2)*exp(-1*(x-m)^2/(2*s2))
}
nxx<-seq(-10, 10, by=0.1)
ny<-normal_Pdf(nxx, 1, 2*2)
nx<-rnorm(1000000, 1,2)
plot(nxx, ny, col=2, lty=2)
lines(density(nx, n=16), col=3)


x<-hist(yy, breaks=2000000,plot=F)#, #prob=TRUE, xlim=c(0,20))
hwidth<-x$breaks[2]-x$breaks[1]
histArea<-hwidth*sum(x$counts)

barHeight<-x$counts/histArea
lines(x$breaks[-length(x$breaks)], barHeight, col=2, lwd=2)

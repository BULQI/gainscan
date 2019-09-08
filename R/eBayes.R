#EBayes module, used to take care the calculation of priors of the variance
#see Smyth 2004 for reference

#
##now start fitting the model to estimate the EBayes parameters
#defining a S3 function for calculating the prior variance and df 
#based on the observed variance

#Oxygen comments
#'@title S3 function to calculate the priors for sample variance
#'@description the function to calculate the priors for sample variance 
#'	following work in Smyth 2004 for reference
#'@details It assume a scaled chisquare prior distribution with d0 and s0
#'	1/si^2~1/(d0*s0^2) chi^2(d0)
#'	check details in Smyth 2004
#'
#'@param x numeric vector of sample variance 
#'@param df numeric vector of degree of freedom for sample standard deviations.
#' 		it is possible to have missing values, where df's for 
#'		different groups vary
#'
#'@return list with d0 (prior df) and s02 (prior variance)

#'@examples
#' 
#' nGene<-1000
#' nTreatment<-2   #number of different beta
#' sampleSize<-100 ##this is the number of repeats for each group

#' alpha.mean<-0  #variance for alpha prior
#' alpha.sigma<-3
#' beta.mean<-0
#' beta.sigma<-2   #variance for beta prior
#'
#' gammaN0.sigma<-10  # the unscaled factor for gamma given gamma<=>0
#' p_diff<-0.01

#' #priors for variance distribution
#' d0<-5
#' s02<-2

#' #call it
#' set.seed(2004);
#' dataExpList<-simulateExpression(nGene, nTreatment, sampleSize,
#'					alpha.mean=alpha.mean, alpha.sigma=alpha.sigma,
#'					beta.mean=beta.mean, beta.sigma=beta.sigma,
#'					prob.nonZero=0.01, gamma.sigma=gammaN0.sigma,
#'					epsilon.d0=d0, epsilon.s0=s0
#'					)
#' #calculate the sample variance from the data
#'	dataExp<-dataExpList[[1]]
#'  
#' #reformat it into a vector
#' s2<-as.vector(sampleVariance(dataExp,nTreatment, sampleSize))
#' df<-as.vector(sampleSizes(dataExp,nTreatment, sampleSize))-1
#' r2<-sampleVariancePrior(s2,df)
#'@export
sampleVariancePrior<-function (x,df)
{
	eg<-log(x)-digamma(df/2)+log(df/2)
	eg_bar<-mean(eg)
	tg_d02<-mean((eg-eg_bar)*(eg-eg_bar)*length(x)/(length(x)-1)-trigamma(df/2))
	
	#here we need to be careful, as was discussed in my notes
	#it could be zeroin the case where we draw from 
	#a constant variance distribution. In this case, it could be negative
	#due to numerical error. So we need to handle this decently
	if(tg_d02 < -1E-0)
	{
		stop("ERROR: negative values for trigammaInver input. Stop!")
	}
	if(tg_d02<0) #in this case, must be larger than -1E-6, treat this as a numerical error
	{
		tg_d02<-0
	}
	
	d0<-trigammaInverse(tg_d02)*2
	d0_num<-d0
	if(is.infinite(d0))
	{
		d0_num<- .Machine$double.xmax
	}
	s02<-exp(eg_bar+digamma(d0_num/2)-log(d0_num/2))
	#inf how ??
	prior<-list(d0=d0, s02=s02)
	
} 
#the function to calculate the sample variance based on data
#Oxygen2 comments
#'@title S3 function to calculate the sample variances based on data
#'@description the function to calculate the sample variances based on the observed data
#'@details It assumes the input data matrix with a format of genes by groups.
#'	It calculates the sample variance for each replicated group. It allows the 
#'	missing values.
#'@param x numeric matrix holding the observed data assuming gene x treatment groups
#'	it also assumes G1R1, G1R2, ..,G1Rn, G2R1,...G2Rn......GmR1...GmRn
#'@param group numeric number of treatment groups. 
#'@param sample.size numeric number of repeats in each group. assuming a blanced design
#'	 for this current version.
#'
#'@return the data matrix contain s2, the sample variance for each groups; 
##  '	mean, the sample mean for each groups
##  ' sample size, the actual sample size by removing the missing values if there are.

#'@examples
#' 
#' nGene<-1000
#' nTreatment<-2   #number of different beta
#' sampleSize<-100 ##this is the number of repeats for each group

#' alpha.mean<-0  #variance for alpha prior
#' alpha.sigma<-3
#' beta.mean<-0
#' beta.sigma<-2   #variance for beta prior
#'
#' gammaN0.sigma<-10  # the unscaled factor for gamma given gamma<=>0
#' p_diff<-0.01

#' #priors for variance distribution
#' d0<-5
#' s0<-2

#' #call it
#' set.seed(2004);
#' dataExpList<-simulateExpression(nGene, nTreatment, sampleSize,
#'					alpha.mean=alpha.mean, alpha.sigma=alpha.sigma,
#'					beta.mean=beta.mean, beta.sigma=beta.sigma,
#'					prob.nonZero=0.01, gamma.sigma=gammaN0.sigma,
#'					epsilon.d0=d0, epsilon.s0=s0
#'					)
#' dataExp<-dataExpList[[1]]
#' #calculate the sample variance from the data
#' svList<-sampleVariance(dataExp,nTreatment, sampleSize)
#' #reformat it into a vector
#' s2<-as.vector(svList)
#' #visualizing the sample variance
#' plot(density(d0*s0^2/rchisq(100000, d0)),col=2, lty=3) 
#' #plot the true distribution 
#' lines(density(s2))
#'@export
sampleVariance<-function(x,group,sample.size)
{
	#first, let's check for the data integrity 
	if(missing(x))
	{
		stop("x, the input data matrix is not correctly set");
	}
	if(class(x)!="matrix")
	{
		stop("x, the input data is not in a correct format");
	}
	
	if(missing(group))
	{
		stop("group, the number of groups is not correctly set");
	}
	if(missing(sample.size))
	{
		stop("sample.size is not correctly set");
	}
	#check for the groups and sample.size
	if(length(sample.size)==1)
	{
		sample.size<-rep(sample.size, group)
	}
	groupInfo<-rep(seq(1:group),times=sample.size)
	tY<-data.frame(t(x))
	tY$group<-groupInfo
	varY<-aggregate(tY,by=list(tY$group),FUN=var, na.rm=TRUE) #by doing aggregate, the new dataframe is arranged to have first column for "by" criteria.
		#that is why in the next step, we need to get rid of the first column and last one, since the last one is the one we added as the group info.
	s2<-varY[,c(-1,-1*(length(x[,1])+2))]
	#missY<-is.na(tY[,c(-1*length(tY[1,]))])
	#missY_dt<-data.frame(1-missY)
	#missY_dt$group<-groupInfo
	#dfY<-aggregate(missY_dt,by=list(missY_dt$group),FUN=sum, na.rm=TRUE)
	#df<-dfY[,c(-1,-1*(length(x[,1])+2))]-1 #minus one since we are one df is used to calculated mean
	t(s2)
	#s2<-as.vector(as.matrix(varK))

	#df<-as.vector(as.matrix(dfK))
	
	#l<-list(s2=s2, df=df)
}

#'@describeIn sampleVariance sampleMean calculate the sample mean based on data
#'
#'@export
sampleMean<-function(x,group,sample.size)
{
	#first, let's check for the data integrity 
	if(missing(x))
	{
		stop("x, the input data matrix is not correctly set");
	}
	if(class(x)!="matrix")
	{
		stop("x, the input data is not in a correct format");
	}
	
	if(missing(group))
	{
		stop("group, the number of groups is not correctly set");
	}
	if(missing(sample.size))
	{
		stop("sample.size is not correctly set");
	}
	#check for the groups and sample.size
	if(length(sample.size)==1)
	{
		sample.size<-rep(sample.size, group)
	}
	groupInfo<-rep(seq(1:group),times=sample.size)
	tY<-data.frame(t(x))
	tY$group<-groupInfo
	meanY<-aggregate(tY,by=list(tY$group),FUN=mean, na.rm=TRUE) #by doing aggregate, the new dataframe is arranged to have first column for "by" criteria.
		#that is why in the next step, we need to get rid of the first column and last one, since the last one is the one we added as the group info.
	m<-meanY[,c(-1,-1*(length(x[,1])+2))]
	#missY<-is.na(tY[,c(-1*length(tY[1,]))])
	#missY_dt<-data.frame(1-missY)
	#missY_dt$group<-groupInfo
	#dfY<-aggregate(missY_dt,by=list(missY_dt$group),FUN=sum, na.rm=TRUE)
	#df<-dfY[,c(-1,-1*(length(x[,1])+2))]-1 #minus one since we are one df is used to calculated mean
	#s2<-as.vector(as.matrix(varK))

	#df<-as.vector(as.matrix(dfK))
	
	#l<-list(s2=s2, df=df)
	t(m)
}

#'@describeIn sampleVariance sampleSizes determine the actual sample size
#'
#'@export
sampleSizes<-function(x,group,sample.size)
{
	#first, let's check for the data integrity 
	if(missing(x))
	{
		stop("x, the input data matrix is not correctly set");
	}
	if(class(x)!="matrix")
	{
		stop("x, the input data is not in a correct format");
	}
	
	if(missing(group))
	{
		stop("group, the number of groups is not correctly set");
	}
	if(missing(sample.size))
	{
		stop("sample.size is not correctly set");
	}
	#check for the groups and sample.size
	if(length(sample.size)==1)
	{
		sample.size<-rep(sample.size, group)
	}
	
	groupInfo<-rep(seq(1:group),times=sample.size)
	tY<-data.frame(t(x))
	tY$group<-groupInfo
	#meanY<-aggregate(tY,by=list(tY$group),FUN=mean, na.rm=TRUE) #by doing aggregate, the new dataframe is arranged to have first column for "by" criteria.
		#that is why in the next step, we need to get rid of the first column and last one, since the last one is the one we added as the group info.
	#m<-meanY[,c(-1,-1*(length(x[,1])+2))]
	missY<-is.na(tY[,c(-1*length(tY[1,]))])
	missY_dt<-data.frame(1-missY)
	missY_dt$group<-groupInfo
	dfY<-aggregate(missY_dt,by=list(missY_dt$group),FUN=sum, na.rm=TRUE)
	df<-dfY[,c(-1,-1*(length(x[,1])+2))] #-1 #<-==here we don't subtract one, since this is the sampleSize, but not the df (minus one since we are one df is used to calculated mean)
	#s2<-as.vector(as.matrix(varK))

	#df<-as.vector(as.matrix(dfK))
	
	#l<-list(s2=s2, df=df)
	t(df)
}

##calculate the posterior sample variance based on the prior and observated sampler variance
#'@title calculate the posterior sample variance mean based data and prior
#'@description by assuming a scaled Chi-square distribution of the sample variances,
#'	it calculated the posterior mean of the variances based on 
#'	the prior and data (observed variances).
#'@details It calculates the posterior mean of the variances based on
#'	the work by Smyth 2004 as below,
#'	s0^2_tilde=(d0*s0^2+d*s^2)/(d0+d)
#'	where s0^2_tilde is the posterior mean of sample variance, 
#'	d0 and s0^2 are the parameters of the prior scaled chi-square distribution
#'	d and s^2 are the observed degree freedom and sample variance.
#'
#'@param d numeric the degree of freedom of the prior distribution
#'@param s2 numeric the variance of the prior distribution
#'@param d0 numeric or vector the degree of freedom of sample data.
#'	It should be one less than the number of data points 
#'@param s02 numeric the sample variance based on the data
#'@return a vector or a numeric as the posterior mean of the sample variance 
#'@seealso sampleVariance sampleVariancePrior rScaledChisq
#'@export
sampleVariancePosterior<-function(d, s2, d0, s02 )
{
	x<-(s2*d+d0*s02)/(d0+d)
	return(x)
}

#transform the data to stabilize the sample variance in order to run linear regression
#'@title transform the data to stabilize the sample variance
#'@description transform the data to keep the mean/location unchanged, but 
#'	stablize the variance.
#'@details to transform the data, we do the following
#'		x'=(x-x_bar)/s+x_bar
#'	where x' is the transformed data, x_bar is the group mean
#'	and s is the posterior mean of sample variance. \code{\link{sampleVariancePosterior}}
#'	First, the input data we assume a data matrix. In order to do the
#'	transformation, we first calcuate the sample variance, and
#'	then fit a scaled chi-square distribution. Finally, with the 
#'	super parameters fitted, we calculated the posterior mean
#' 	of sample variance of each group. In the last step, we
#'	transform the data
#'@param x matrix contains the data in a form of protein by group. Each
#'	row of the data is one replicate of one group. Also all the replicates
#'	are grouped together adjacently. 
#'@param group numeric the number of groups for the data
#'@param sample.size numeric the number of replicates in each group
#'@return numeric matrix with the transformed data which has been stablized
#'	for its variances
#'@export
transformData<-function(x,group, sample.size)
{
	#first, let's check for the data integrity 
	if(missing(x))
	{
		stop("x, the input data matrix is not correctly set");
	}
	if(class(x)!="matrix")
	{
		stop("x, the input data is not in a correct format");
	}
	
	if(missing(group))
	{
		stop("group, the number of groups is not correctly set");
	}
	if(missing(sample.size))
	{
		stop("sample.size is not correctly set");
	}
	#check for the groups and sample.size
	if(length(sample.size)==1)
	{
		sample.size<-rep(sample.size, group)
	}
	#now we should have sample.size as a vector
	ncolum<-dim(x)[2]
	
	if(ncolum!=sum(sample.size))
	{
		stop("the number of columns doesn't equal to the sample description by the input, please check!")
	}
	#now everything looks is correct
	
	#get the sample mean and variance
	sMean<-sampleMean(x, group, sample.size)
	sVar<-sampleVariance(x,group,sample.size)
	sdf<-sampleSizes(x,nTreatment, sampleSize)-1
	#get the sample variance prior by assuming a scale-chisquare
	s2<-as.vector(as.matrix(sVar))
	df<-as.vector(as.matrix(sdf))
	
	r2<-sampleVariancePrior(s2,df)
	s02<-r2[[2]]
	df0<-r2[[1]]
	
	#calculate the posterior variance
	var_pos<-sampleVariancePosterior(sdf,sVar,  df0, s02 )
	
	#prepare the mean for ".matrix" operation
	#mean_pos<-x
	
	#now do the transformation
	running_index<-0
	for(i in c(1:group))
	{
		#cat("i:",i, "\n")
		for(j in c(1:sample.size[i]))
		{
		
			#cat("\tj:",j, "\n")
			#mean_pos[,running_index]<-sMean[,i]
				x[,running_index+j]<-(x[,running_index+j]-sMean[,i])/sqrt(var_pos[,i])+sMean[,i]
		}
		running_index<-running_index+sample.size[i];
	}
	#now do the transformation with ".matrix" operation
	#x<-(x-sMean)/sqrt(var_pos)+sMean;
	return(x)
}

#the code used to generate a sample (variance) following a scaled chi-Square d'n
#
##oxygen2 comments
#'@title Scaled Chi-Square distribution
#'@description generate samples from the Scaled Chi-Square distrubtion.
#'  This distribution is used to describe the sample
#'	variances with parameters d0 and s02.
#'@details This is described in details in Smyth 2004 work. It has the following
#'	distribution
#'	\deqn{d0*s0^2/s^2~chi^2(d0)}
#' or
#'	\deqn{1/s^2~d0*s0^2*chi^2(d0)}
#'
#'@param n integer the number of samples to generated
#'@param df0 interger the degree of freedom of the distribution
#'@param s02 numeric the prior sample variance.
#'
#'@return a vector of sample variances
#'
#'@export 
rScaledChisq<-function(n, df0, s02)
{
	if(missing(n))
	{
		stop("ERROR: argument 'n' is missing, with no default. ");
	}
	if(missing(df0))
	{
		stop("ERROR: argument 'df' is missing, with no default. ");
	}
	if(missing(s02))
	{
		stop("ERROR: argument 's0' is missing, with no default. ");
	}
	#now do it.
	if(n<=0)
	{
		stop("ERROR: invalid argument.");
	}
	if(df0<=0)
	{
		stop("ERROR: invalid argument.");
	}
	if(s02<=0)
	{
		stop("ERROR: invalid argument.");
	}
	#first call the chisq
	x<-rchisq(n,df0);
	x<-df0*s02/x;
	return(x)
}
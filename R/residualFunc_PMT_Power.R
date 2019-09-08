###residualFunc_PMT_Power. 
###In this file, we define all the residual functions for PMT Power model fitting
####for gainAdjust.R project.
###started 1/27/2019 Feng @ Boston University

pmt_saturated<-2^16-1;
pmt_baseline<-log(0.1)

##we are looking for a solution to overcome the issue where we run out of limit

	#'@title the PMT power function with baseline
	#'@description This is the function to describe the PMT out intensity based on the
	#'gain and input light flux. It is a power function of the dynode voltage with 
	#'a baseline. 
	#'@details The function has the following format \cr
	#' log(y)=log(exp(B)+exp(delta*(log(x)+phi)))
	#'\cr
	#'It has a maximum value. On a 16-bit representation, it is 2^16-1.
	#'The function takes in the parameter and the input voltages. 
	#'The the function values is the final fluoreoutput. In this implementation,
	#'phi is assumed to be log-transformed. This way, we enforce phi is a positive
	#'number.
	#'@param pars a vector of 4 numerics, B, G, phi and delta. phi and B are log-transformed 
	#'to force them to be non-negative in the model. exp(B) >0 exp(phi)>0 
	#'@param x a vector of numerics. 
	#'@param bool a string to indicate whether the saturatoin will be done. TRUE by default.
	#'@return  a vector of numerics of similar length as x. 
	#'@export
	
	#this one is for single parameters but a vector of x/voltage
	f_PowerBaselinePMT<-function(pars,x, saturation=TRUE)
	{
		if(length(pars)!=4){
			stop("pars are not correctly specified!")
		}
		B<-pars[1]
		G<-pars[2]
		phi<-pars[3]
		delta<-pars[4]
		#g<-pars[5]
		
		#y<-x; #make the holder first
		#y[]<- -1
		##split the array into to parts
		#cat("B:",B, ";delta:",delta, "phi:",phi,"===");
		#y<-log(B+exp(delta*log(x*phi)))  #first try, phi could be negative. that doesn't make sense
		#y<-log(B+exp(delta*(log(x)+phi))) #secont try, phi is log(phi') now, force it to be positive.
		y<-f_PowerBaselinePMT_core(B, G,phi,delta,x, saturation);
		#log(exp(B)+exp(delta*(log(x)+phi))) #now B=log(B), for exp(B) to be 
		#cat("y:",y,"-----x:",x);
		#cat("\n");
		#y<-exp(y);
		#if(saturation){
		#	y[y>G]<-G;
		#}
		return(y)
	}
	
	#this is the function version that takes 
	#'@export
	f_PowerBaselinePMT_phi<-function(B, G,phi,delta,x, saturation=TRUE)
	{
		if(missing(B))
		{
			stop("please specify B....!!!")
		}
		if(missing(G))
		{
			stop("please specify G....!!!")
		}
		if(missing(phi))
		{
			stop("please specify phi....!!!")
		}
		if(missing(delta))
		{
			stop("please specify delta....!!!")
		}
		if(missing(x))
		{
			stop("please specify x....!!!")
		}
		if(length(B)>1)
		{
			B<-B[1]
			warning("more than one elements are specified for B. Only the first one is used!!")
		}
		if(length(G)>1)
		{
			G<-G[1]
			warning("more than one elements are specified for G. Only the first one is used!!")
		}
		if(length(delta)>1)
		{
			delta<-delta[1]
			warning("more than one elements are specified for delta. Only the first one is used!!")
		}
		if(length(x)>1)
		{
			x<-x[1]
			warning("more than one elements are specified for x. Only the first one is used!!")
		}
		y<-f_PowerBaselinePMT_core(B, G,phi,delta,x, saturation);
		#log(exp(B)+exp(delta*(log(x)+phi))) #now B=log(B), for exp(B) to be 
		#cat("y:",y,"-----x:",x);
		#cat("\n");
		#y<-exp(y);
		#if(saturation){
		#	y[y>G]<-G;
		#}
		return(y)
	
	}
	#this is the function version that takes 
	#'@export
	#the caller will make sure the parameters are correctly set up:
	#importantly, one of phi or x could be a vector at one time. others can 
	# only be single parameter.
	f_PowerBaselinePMT_core<-function(B, G,phi,delta,x, saturation=TRUE)
	{
		#if(length(pars)!=){
		#	stop("pars are not correctly specified!")
		#}
		#B<-pars[1]
		#G<-pars[2]
		#phi<-pars[3]
		#delta<-pars[4]
		#g<-pars[5]
		if(missing(B))
		{
			stop("please specify B....!!!")
		}
		if(missing(G))
		{
			stop("please specify G....!!!")
		}
		if(missing(phi))
		{
			stop("please specify phi....!!!")
		}
		if(missing(delta))
		{
			stop("please specify delta....!!!")
		}
		if(missing(x))
		{
			stop("please specify x....!!!")
		}
		
		#y<-phi; #make the holder first
		#y[]<- -1
		##split the array into to parts
		#cat("B:",B, ";delta:",delta, "phi:",phi,"===");
		#y<-log(B+exp(delta*log(x*phi)))  #first try, phi could be negative. that doesn't make sense
		#y<-log(B+exp(delta*(log(x)+phi))) #secont try, phi is log(phi') now, force it to be positive.
		y<-log(exp(B)+exp(delta*(log(x)+phi))) #now B=log(B), for exp(B) to be 
		#cat("y:",y,"-----x:",x);
		#cat("\n");
		y<-exp(y);
		if(saturation){
			y[y>G]<-G;
		}
		return(y)
	}


####always assume the first set of x and y are the reference ones
####meaning without be adjusted for shift parameter k
#'@title residual function for the PMT power model with baseline 
#'@description this is the residual function for fitting the PMT power model with baseline. It 
#'	assumes a log residual error, but could be switched to do non-log residuals.
#'@param pars a vector holding the parameters for the model. we assume the following order
#'	c(B, delta, phi's). The total length of pars is length of(x)+2. phi and B is assumed
#' to be log-transformed.
#'@param y data matrix holding the intensity. Arranged as gene by voltage.
#'@param x is the vector holding the input voltages. It has the identical length as 
#'	for the columns of y.
#'@param G the maximum possible intensity. In a 16 bit presentaion, it has
#'	a value of 2^16-1.
#'@param residual.mode indicating whether to use a log or regular scale of residuals.
#'@return a vector of the residual. It has a total length of number of y elements.
#'
#'@export
gainAdjust_fnResPbp<-function(pars, #phi has to be log transformed.
		y, #the response, assuming log transformed
		x, #the independent variable, log-transformed
		G=pmt_saturated,
		residual.mode=c("log","regular")
		#xmid, #x value at the middel point of y
		#scal, #slope at xmid
		#g #the asymetric factor
		)
	{
		#cat("residua.........fun");
		if(class(y)!="matrix")
		{
			stop("input y is not specified as matrix, please check");
		}
		if(class(x)!="matrix"&&class(x)!="numeric")
		{
			stop("input x is not specified as matrix or numeric, please check");
		}
		
		if(class(x)=="numeric")
		{
			if(dim(y)[2]!=length(x))
			{
				stop("x length doesn't match y dimension. please check");
			}
			
		}
		residual.mode<-match.arg(residual.mode);
		#now we need to get the parameters out from pars
		B<-pars[1];
		delta<-pars[2];
		phi<-pars[-c(1,2)];
		
		#check to make sure, we have enough phi'scal
		if(length(phi)!=dim(y)[1])
		{
			stop("the input parameters are not enought for fitting the data"); 
		}
		####check the input to make sure the input are in the correct format
		res<-rep(0, dim(y)[1]*dim(y)[2])
		
		for( i in 1:dim(y)[1])
		{
			switch(residual.mode,
				log={
					res[dim(y)[2]*(i-1)+c(1:dim(y)[2])]<-(log(y[i,])-log(f_PowerBaselinePMT(c(B, G, phi[i], delta),x)));#/log(x);
					},
				regular={
					res[dim(y)[2]*(i-1)+c(1:dim(y)[2])]<-y[i,]-f_PowerBaselinePMT(c(B, G, phi[i], delta),x);}
			)
		}
		
		#now adding variance
		#res<-res; #res/x
		#pred<-f5pl(c(a,d,xmid, scal,g),xinput);
		return(res);
	}
	
	#'@title residual function for the PMT power model with baseline 
#'@description this is the residual function for fitting the PMT power model with baseline
#'	for SINGLE gene data series. It 
#'	assumes a log residual error, but could be switched to do non-log residuals.
#'@details Here we based on the fixed B, G and delta. We fit the phi parameter for one 
#	gene data series. 
#'@param pars a numeric holding the parameter phi for the model. we assume the following order
#'	phi is assumed to be log-transformed.
#'@param y data matrix holding the intensities for one gene with varying PMT voltages. 
#'	The length of it is identical to the x vector
#'@param x is the vector holding the input voltages. It has the identical length as 
#'	for the columns of y.
#'@param G the maximum possible intensity. In a 16 bit presentaion, it has
#'	a value of 2^16-1.
#'@param B the baseline level, assumed to be log-transformed in the model.
#'@param delta the exponent of the PMT power model
#'@param residual.mode indicating whether to use a log or regular scale of residuals.
#'@return a vector of the residual. It has a total length of number of y elements.
#'
#'@export
gainAdjust_fnResPbp_phi<-function(pars, 
		y, #the response, assuming log transformed
		x, #the independent variable, log-transformed
		B, G=pmt_saturated, delta,
		residual.mode=c("log","regular")
		#xmid, #x value at the middel point of y
		#scal, #slope at xmid
		#g #the asymetric factor
		)
	{
		#cat("residua.........fun");
		if(class(y)!="numeric")
		{
			stop("input y is not specified as matrix, please check");
		}
		
		if(class(x)=="numeric")
		{
			if(length(y)!=length(x))
			{
				stop("x length doesn't match y dimension. please check");
			}
			
		}
		if(missing(B))
		{
			stop("INput B is not specified, please check!!");
		}
		if(missing(delta))
		{
			stop("input delta is not specified, please check!!!");
		}
		residual.mode<-match.arg(residual.mode);
		#now we need to get the parameters out from pars
		#B<-pars[1];
		#delta<-pars[2];
		phi<-pars; #[-c(1,2)];
		
		#check to make sure, we have enough phi'scal
		#if(length(phi)!=length(y))
		#{
		#	stop("the input parameters are not enought for fitting the data"); 
		#}
		####check the input to make sure the input are in the correct format
		#res<-rep(0, length(y))
		
		#for( i in 1:dim(y)[1])
		#{
			switch(residual.mode,
				log={
					res<-(log(y)-log(f_PowerBaselinePMT(c(B, G, phi, delta),x)));#/log(x);#(log(f_PowerBaselinePMT(c(B, G, phi, delta),x))+5);
					},
				regular={
					res<-y-f_PowerBaselinePMT(c(B, G, phi, delta),x);}
			)
		#}
		#pred<-f5pl(c(a,d,xmid, scal,g),xinput);
		return(res);
	}#end of 
	
#--->gainAdjust.fnRes5pFpl<-function( pars, #parameters 

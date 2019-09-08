###residualFunc.R
###In this file, we define all the residual functions for 5PL fitting
####for gainAdjust.R project.
###started 8/7/2017 Feng @ Boston University


##--->	Don't export
####function NOT working, need manually change to work!!!
# it assumes 
#===>	
##we are looking for a solution to overcome the issue where we run out of limit
## when calculate exp((xmid-x)/scal)
## the solution is like this,
## take a exp(log((d-a)/(1+exp((xmid-x)/scal))^g))
##   rearrange this  exp(log(d-a)-log((1+exp(xmid-x)/scal)^g))
##					exp(log(d-a)-g*log(1+exp(xmid-x)/scal))
##		there are one case we need to consider, when exp is too big or exp is too small (both relative to 1).
##		where exp() is too big or (xmid-x)>709, we approximate to get log(1+exp(xmid-x)/scal)~=(xmid-x)/scal 
##		the whole thing becomes exp(log(d-a)-g*(xmid-x)/scal)
	#the five parameter logistic function
	#take in the parameter 5 parameters and get the function values
	#'@export
	f5pl<-function(pars,x)
	{
		if(length(pars)!=5){
			stop("pars are not correctly specified!")
		}
		a<-pars[1]
		d<-pars[2]
		xmid<-pars[3]
		scal<-pars[4]
		g<-pars[5]
		
		y<-x; #make the holder first
		y[]<- -1
		#split the array into to parts
		regular<-which((xmid-x)/scal<=709)
		if(length(regular)==0){
			y<-a+exp(log(d-a)-g*(xmid-x)/scal)
		} else {
			y[regular]<-a+(d-a)/(1+exp((xmid-x[regular])/scal))^g
			y[-regular]<-a+exp(log(d-a)-g*(xmid-x[-regular])/scal)	
		}
		return(y)
		#pred<-a+(d-a)/(1+exp((xmid-x)/scal))^g
	}


####always assume the first set of x and y are the reference ones
####meaning without be adjusted for shift parameter k
gainAdjust.fplFit<-function(
		y, #the response, assuming log transformed
		x, #the independent variable, log-transformed
		k, #vector of k to changing x 
		k2,
		#a, #par a, smallest y log
		#d, #par d, largest y log
		xmid, #x value at the middel point of y
		scal, #slope at xmid
		g #the asymetric factor
		)
	{
		####check the input to make sure the input are in the correct format
		if(length(x)!=length(y))
		{
			stop("input x and y are not equal in length")
		}
		xinput<-x+c(rep(0,length(x)/3),rep(k,length(x)/3),rep(k2,length(x)/3))
		
		#xinput<-x+c(rep(0,length(x)/3),rep(k[1],length(x)/3),rep(k[2],length(x)/3))
		#pred<-a+(d-a)/(1+exp((xmid-xinput)/scal))^g   #<----changed this
		pred<-f5pl(c(a,d,xmid, scal,g),xinput);
		(y-pred)
	}
	
	##==>
	##careful, this function ASSUMES we HAVE SHIFTED the x input!!!<--------
	##-->
	#pars=c(a,d, xmid, scal, g), no other parameters
	#aggregate
	#'@export
gainAdjust.fnRes5pFpl<-function( pars, #parameters 
		y, #the response, assuming log transformed
		x, #the independent variable, log-transformed
		#k, #vector of k to changing x 
		#a, #par a, smallest y log
		#d, #par d, largest y log
		#xmid, #x value at the middel point of y
		#scal, #slope at xmid
		ylog=TRUE,
		#aggregated=TRUE #indicated whether the data has been aggregated, mean the two repeats has been averaged
	#				by default it is not. So in this case we need to give the two repeat the identical shift, k
		model.weight="log", order=1 #order is used for generation of exponential weight matrix, the order of the power. Not used if other type of weight matrix is applied.
		)
	{
		####check the input to make sure the input are in the correct format
		if(length(x)!=length(y))
		{
			stop("input x and y are not equal in length")
		}
		
		#now we need to figure out the correct xinput
		#the pars is a vector and the first 3 always xmid scal, and g.
		#the rest are ks
		#ks<-pars[4:length(pars)]
		a<-pars[1]
		d<-pars[2]
		xinput<-x
		xmid<-pars[3]
		scal<-pars[4]
		g<-pars[5]
		#if(ylog)
		#{
		#	a<-log(a)
		#	d<-log(d)
		#}
		#pred<-a+(d-a)/(1+exp((xmid-xinput)/scal))^g
		pred<-f5pl(c(a,d,xmid, scal,g),xinput);
		m<-weight.matrix(pred,model.weight,order)
		(y-pred)/m  #log(abs(pred))#*pred#^2#*(pred)
	}	
	
	##==>
	##careful, this function ASSUMES we HAVE SHIFTED the x input!!!<--------
	##-->
	#pars=c(a,d, xmid, scal, g), no other parameters
	#aggregate
	#'@export
gainAdjust.fnRes4pFpl<-function( pars, #parameters 
		y, #the response, assuming log transformed
		x, #the independent variable, log-transformed
		#k, #vector of k to changing x 
		#a, #par a, smallest y log
		d, #par d, largest y log
		#xmid, #x value at the middel point of y
		#scal, #slope at xmid
		ylog=TRUE,
		#aggregated=TRUE #indicated whether the data has been aggregated, mean the two repeats has been averaged
	#				by default it is not. So in this case we need to give the two repeat the identical shift, k
		model.weight="sqrt", order=1 #order is used for generation of exponential weight matrix, the order of the power. Not used if other type of weight matrix is applied.
		)
	{
		####check the input to make sure the input are in the correct format
		if(length(x)!=length(y))
		{
			stop("input x and y are not equal in length")
		}
		
		#now we need to figure out the correct xinput
		#the pars is a vector and the first 3 always xmid scal, and g.
		#the rest are ks
		#ks<-pars[4:length(pars)]
		#a<-pars[1]
		#d<-pars[2]
		xinput<-x
		a<-pars[1]
		xmid<-pars[2]
		scal<-pars[3]
		g<-pars[4]
		if(missing(d))
		{
			d<-pmt_saturated_reading
			if(ylog)
			{
				#a<-log(a)
				d<-log(d)
			}
		}
		#
		#pred<-a+(d-a)/(1+exp((xmid-xinput)/scal))^g
		pred<-f5pl(c(a,d,xmid, scal,g),xinput);
		m<-weight.matrix(pred,model.weight,order)
		(y-pred)/m  #sqrt(abs(pred))#/log(abs(pred))#*pred#^2#*(pred)
		#(y-pred)
	}	
	##==>
	##careful, this function ASSUMES we HAVE SHIFTED the x input!!!<--------
	##-->
	#pars=c(xmid, scal, g), no other parameters
	#aggregate
	#'@export
gainAdjust.fnRes3pFpl<-function( pars, #parameters 
		y, #the response, assuming log transformed
		x, #the independent variable, log-transformed
		#k, #vector of k to changing x 
		a, #par a, smallest y log
		d, #par d, largest y log
		#xmid, #x value at the middel point of y
		#scal, #slope at xmid
		ylog=TRUE,
		#aggregated=TRUE #indicated whether the data has been aggregated, mean the two repeats has been averaged
	#				by default it is not. So in this case we need to give the two repeat the identical shift, k
		model.weight="sqrt", order=1 #order is used for generation of exponential weight matrix, the order of the power. Not used if other type of weight matrix is applied.
		)
	{
		####check the input to make sure the input are in the correct format
		if(length(x)!=length(y))
		{
			stop("input x and y are not equal in length")
		}
		
		#now we need to figure out the correct xinput
		#the pars is a vector and the first 3 always xmid scal, and g.
		#the rest are ks
		#ks<-pars[4:length(pars)]
		#a<-pars[1]
		#d<-pars[2]
		xinput<-x
		#a<-pars[1]
		xmid<-pars[1]
		scal<-pars[2]
		g<-pars[3]
		if(missing(a))
		{
			a<-pmt_lowest_reading
			if(ylog)
			{
				a<-log(a)
			}
		}
		if(missing(d))
		{
			d<-pmt_saturated_reading
			if(ylog)
			{
				a<-log(a)
			}
		}
		#if(ylog)
		#{
		#	#a<-log(a)
		#	d<-log(d)
		#}
		#pred<-a+(d-a)/(1+exp((xmid-xinput)/scal))^g
		pred<-f5pl(c(a,d,xmid, scal,g),xinput);
		m<-weight.matrix(pred,model.weight,order)
		(y-pred)/m  #sqrt(abs(pred))#/log(abs(pred))#*pred#^2#*(pred)
		#(y-pred)
	}	
	
	##Residual function used by nlsLM it will estimate the
	## ks for shifting/aligning x input.
#pars=c(xmid, scal, g, k1, k2,...)
#aggregate
#'@export
gainAdjust.fnRes3pShift<-function( pars, #parameters 
		y, #the response, assuming log transformed
		x, #the independent variable, log-transformed
		#k, #vector of k to changing x 
		a, #par a, smallest y log
		d, #par d, largest y log
		#xmid, #x value at the middel point of y
		#scal, #slope at xmid
		#g #the asymetric factor
		xlen, #the length of distinct xs, x<-c(250,300,350,400,450,500,550,600,650), then xlens=9
		shift.index=1, #this is the index of the data set in the y array, to which all other data series shift.
		ylog=TRUE,
		aggregated=TRUE #indicated whether the data has been aggregated, mean the two repeats has been averaged
	#				by default it is not. So in this case we need to give the two repeat the identical shift, k
		, model.weight="sqrt", order=1 #order is used for generation of exponential weight matrix, the order of the power. Not used if other type of weight matrix is applied.
		)
	{
		####check the input to make sure the input are in the correct format
		if(length(x)!=length(y))
		{
			stop("input x and y are not equal in length")
		}
		
		#now we need to figure out the correct xinput
		#the pars is a vector and the first 3 always xmid scal, and g.
		#the rest are ks
		ks<-pars[4:length(pars)]
		if(aggregated)
		{
			if(length(x)/xlen!=length(ks)+1)
				{
					stop("the parameters k for shifting the data x series are not set with correct size!")
				}
		}
		else
		{
			if(length(x)/2/xlen!=length(ks)+1)
				{
					stop("the parameters k for shifting the data x series are not set with correct size--error02!")
				}
		}
		
		if(aggregated)
		{
			xRepeat<-1;
		}
		else
		{
			xRepeat<-2
		}
		#ks<-c(1,2,3)
		#ks<-gainAdjust.insertAt(ks,0,shift.index);
		##--------------v
		if(shift.index==1)
		{
			ks<-c(0,ks)
		} else if(shift.index==length(ks)+1)
		{
			ks<-c(ks,0)
		}else
		{
			#insert
			ks<-c(ks[1:shift.index-1],0,ks[shift.index:length(ks)])
		}
		###-------------------------^
		#ks<-c(0,ks)
		ks<-rep(ks,rep(xRepeat,length(ks)))
		ks<-rep(ks,rep(xlen, length(ks)))
		xinput<-x+ks
		xmid<-pars[1]
		scal<-pars[2]
		g<-pars[3]
		if(missing(a))
		{
			a<-pmt_lowest_reading
			if(ylog)
			{
				a<-log(a)
				#d<-log(d)
			}
		}
		if(missing(d))
		{
			d<-pmt_saturated_reading
			if(ylog)
			{
				#a<-log(a)
				d<-log(d)
			}
		}
		
		#pred<-a+(d-a)/(1+exp((xmid-xinput)/scal))^g
		pred<-f5pl(c(a,d,xmid, scal,g),xinput);
		m<-weight.matrix(pred,model.weight,order)
		(y-pred)/m #sqrt(abs(pred))#^2#*(pred)
	}
	
	
	##Residual function used by nlsLM it will estimate the
	## ks for shifting/aligning x input.
	#IMPORTANTLY, in this 4pL, we allow lowest value param, a, to vary, keeping highest value param,d, fixed 
#pars=c(xmid, scal, g, k1, k2,...)
#aggregate
#'@export
##
gainAdjust.fnRes4pShift<-function( pars, #parameters 
		y, #the response, assuming log transformed
		x, #the independent variable, log-transformed
		#k, #vector of k to changing x 
		#a, #par a, smallest y unlogged
		d, #par d, largest y unlogged
		#xmid, #x value at the middel point of y
		#scal, #slope at xmid
		#g #the asymetric factor
		xlen, #the length of distinct xs, x<-c(250,300,350,400,450,500,550,600,650), then xlens=9
		ylog=TRUE,
		aggregated=TRUE #indicated whether the data has been aggregated, mean the two repeats has been averaged
	#				by default it is not. So in this case we need to give the two repeat the identical shift, k
		, model.weight="sqrt", order=1 #order is used for generation of exponential weight matrix, the order of the power. Not used if other type of weight matrix is applied.
		)
	{
		####check the input to make sure the input are in the correct format
		if(length(x)!=length(y))
		{
			stop("input x and y are not equal in length")
		}
		
		#now we need to figure out the correct xinput
		#the pars is a vector and the first 3 always xmid scal, and g.
		#the rest are ks
		ks<-pars[5:length(pars)]
		if(aggregated)
		{
			if(length(x)/xlen!=length(ks)+1)
				{
					stop("the parameters k for shifting the data x series are not set with correct size!")
				}
		}
		else
		{
			if(length(x)/2/xlen!=length(ks)+1)
				{
					stop("the parameters k for shifting the data x series are not set with correct size--error02!")
				}
		}
		
		if(aggregated)
		{
			xRepeat<-1;
		}
		else
		{
			xRepeat<-2
		}
		#ks<-c(1,2,3)
		ks<-c(0,ks)
		ks<-rep(ks,rep(xRepeat,length(ks)))
		ks<-rep(ks,rep(xlen, length(ks)))
		xinput<-x+ks
		a<-pars[1]
		xmid<-pars[2]
		scal<-pars[3]
		g<-pars[4]
		if(missing(d))
		{
			d<-pmt_saturated_reading
		}
		if(ylog)
		{
			#a<-log(a)
			d<-log(d)
		}
		#pred<-a+(d-a)/(1+exp((xmid-xinput)/scal))^g
		pred<-f5pl(c(a,d,xmid, scal,g),xinput);
		m<-weight.matrix(pred,model.weight,order)
		(y-pred)/m #sqrt(abs(pred))#pred #^2#*(pred)
	}
	
	
	##Residual function used by nlsLM it will estimate the
	## ks for shifting/aligning x input.
	#IMPORTANTLY, in this 4pL, we allow lowest value param, a, to vary, keeping highest value param,d, fixed 
#pars=c(xmid, scal, g, k1, k2,...)
#aggregate
#
#'@export
##
gainAdjust.fnRes5pShift<-function( pars, #parameters 
		y, #the response, assuming log transformed
		x, #the independent variable, log-transformed
		#k, #vector of k to changing x 
		#a, #par a, smallest y unlogged
		#d, #par d, largest y unlogged
		#xmid, #x value at the middel point of y
		#scal, #slope at xmid
		#g #the asymetric factor
		xlen, #the length of distinct xs, x<-c(250,300,350,400,450,500,550,600,650), then xlens=9
		ylog=TRUE,
		aggregated=TRUE #indicated whether the data has been aggregated, mean the two repeats has been averaged
	#				by default it is not. So in this case we need to give the two repeat the identical shift, k
		, model.weight="sqrt", order=1 #order is used for generation of exponential weight matrix, the order of the power. Not used if other type of weight matrix is applied.
		)
	{
		####check the input to make sure the input are in the correct format
		if(length(x)!=length(y))
		{
			stop("input x and y are not equal in length")
		}
		
		#now we need to figure out the correct xinput
		#the pars is a vector and the first 3 always xmid scal, and g.
		#the rest are ks
		ks<-pars[6:length(pars)]
		if(aggregated)
		{
			if(length(x)/xlen!=length(ks)+1)
				{
					stop("the parameters k for shifting the data x series are not set with correct size!")
				}
		}
		else
		{
			if(length(x)/2/xlen!=length(ks)+1)
				{
					stop("the parameters k for shifting the data x series are not set with correct size--error02!")
				}
		}
		
		if(aggregated)
		{
			xRepeat<-1;
		}
		else
		{
			xRepeat<-2
		}
		#ks<-c(1,2,3)
		ks<-c(0,ks)
		ks<-rep(ks,rep(xRepeat,length(ks)))
		ks<-rep(ks,rep(xlen, length(ks)))
		xinput<-x+ks
		a<-pars[1]
		d<-pars[2]
		xmid<-pars[3]
		scal<-pars[4]
		g<-pars[5]
		#if(missing(d))
		#{
		#	d<-pmt_saturated_reading
		#}
		#if(ylog)
		#{
			#a<-log(a)
		#	d<-log(d)
		#}
		#pred<-a+(d-a)/(1+exp((xmid-xinput)/scal))^g
		pred<-f5pl(c(a,d,xmid, scal,g),xinput);
		m<-weight.matrix(pred,model.weight, order)
		(y-pred)/m #sqrt(abs(pred))#pred #^2#*(pred)
	}
	
	##Residual function used by nlsLM it will estimate the
	## ks for shifting/aligning x input.
	#IMPORTANTLY, in this 4pL, we allow lowest value param, a, to vary, keeping highest value param,d, fixed 
#pars=c(xmid, scal, g, k1, k2,...)
#aggregate
#'@export
##
gainAdjust.fnResShift<-function( pars, #parameters 
		y, #the response, 
		x, #the independent variable, log-transformed
		#k, #vector of k to changing x 
		a, #par a, smallest y unlogged
		d, #par d, largest y unlogged
		xmid, #x value at the middel point of y
		scal, #slope at xmid
		g ,#the asymetric factor
		xlen, #the length of distinct xs, x<-c(250,300,350,400,450,500,550,600,650), then xlens=9
		ylog=TRUE, shift.index=1,
		aggregated=TRUE #indicated whether the data has been aggregated, mean the two repeats has been averaged
	#				by default it is not. So in this case we need to give the two repeat the identical shift, k
		, model.weight="sqrt", order=1 #order is used for generation of exponential weight matrix, the order of the power. Not used if other type of weight matrix is applied.
		)
	{
		####check the input to make sure the input are in the correct format
		if(length(x)!=length(y))
		{
			stop("input x and y are not equal in length")
		}
		if(missing(xmid)|missing(scal)|missing(g))
		{
			stop("please specify the parameters (xmid/xscal/g) for 5PL")
		}
		#now we need to figure out the correct xinput
		#the pars is a vector and the first 3 always xmid scal, and g.
		#the rest are ks
		ks<-pars
		if(aggregated)
		{
			if(length(x)/xlen!=length(ks)+1)
				{
					stop("the parameters k for shifting the data x series are not set with correct size!")
				}
		}
		else
		{
			if(length(x)/2/xlen!=length(ks)+1)
				{
					stop("the parameters k for shifting the data x series are not set with correct size--error02!")
				}
		}
		
		if(aggregated)
		{
			xRepeat<-1;
		}
		else
		{
			xRepeat<-2
		}
		#ks<-c(1,2,3)
		#ks<-c(0,ks)
		#ks<-c(1,2,3)
		if(shift.index==1)
		{
			ks<-c(0,ks)
		} else if(shift.index==length(ks)+1)
		{
			ks<-c(ks,0)
		}else
		{
			#insert
			ks<-c(ks[1:shift.index-1],0,ks[shift.index:length(ks)])
		}
		ks<-rep(ks,rep(xRepeat,length(ks)))
		ks<-rep(ks,rep(xlen, length(ks)))
		xinput<-x+ks
		
		if(missing(d))
		{
			d<-pmt_saturated_reading
			if(ylog)
			{
				a<-log(a)
				#d<-log(d)
			}
		}
		if(missing(a))
		{
			a<-pmt_lowest_reading
			if(ylog)
			{
				#a<-log(a)
				d<-log(d)
			}
		}

		
		#pred<-a+(d-a)/(1+exp((xmid-xinput)/scal))^g
		pred<-f5pl(c(a,d,xmid, scal,g),xinput);
		#m<-weight.matrix(pred,model.weight,order)
		#cat("\ny:",y, "\n")
		#cat("\npred:", pred, "\n")
		#cat("ks:", ks, "\n")
		#cat("RSS:", sum((log(y)-log(pred))*(log(y)-log(pred))),"\n")
		log(y)-log(pred) #/ #sqrt(abs(pred))#pred #^2#*(pred)
	}
	
	##Residual function used by nlsLM to fit for individual data set
	## ks for shifting/aligning x input.
	#IMPORTANTLY, in this 5pL, we allow the parametes to be fed in 
#pars=c(k) common k for shifting two duplicated points
#aggregate
##'@export
##

#'@export
disaggregateVec<-function(m, xlen)
{
	if(missing(m)||missing(xlen))
	{
		stop("missing input m/xlen!!")
	}
	if(floor(length(m)/xlen)*xlen!=length(m))
	{
		stop("the length of m is not of xlen fold, please check")
	}
	m.out<-rep(0, 2*length(m))
	n.xlen<-length(m)/xlen
	for(i in 1:n.xlen)
	{
		m.out[c(1:xlen)+(i*2-2)*xlen]<-m[c(1:xlen)+(i-1)*xlen]
		m.out[c(1:xlen)+(i*2-1)*xlen]<-m[c(1:xlen)+(i-1)*xlen]
	}
	m.out
}
#'@export
gainAdjust.fnResShiftSingle<-function( pars, #parameters 
		y, #the response, two duplicated data in vector
		x, #the independent variable, log-transformed
		#k, #vector of k to changing x 
		a, #par a, smallest y unlogged
		d, #par d, largest y unlogged
		xmid, #x value at the middel point of y
		scal, #slope at xmid
		g ,#the asymetric factor
		#xlen, #the length of distinct xs, x<-c(250,300,350,400,450,500,550,600,650), then xlens=9
		#ylog=TRUE, shift.index=1,
		aggregated=TRUE #indicated whether the data has been aggregated, mean the two repeats has been averaged
	#				by default it is not. So in this case we need to give the two repeat the identical shift, k
		, model.weight="uniform", order=1, #order is used for generation of exponential weight matrix, the order of the power. Not used if other type of weight matrix is applied. 
		mvar###this is used to pass in the 
		)
	{
		####check the input to make sure the input are in the correct format
		if(length(x)!=length(y))
		{
			if(length(x)!=length(y)/2)
			{
				stop("input x and y are not equal in length")
			}
			x<-rep(x,2)
			
		}
		if(missing(xmid)|missing(scal)|missing(g))
		{
			stop("please specify the parameters (xmid/xscal/g) for 5PL")
		}
		xlen<-length(x)
		#now we need to figure out the correct xinput
		#the pars is a vector and the first 3 always xmid scal, and g.
		#the rest are ks
		ks<-pars
		#if(aggregated)
		#{
			#if(length(x)/xlen!=length(ks)+1)
			#	{
			#		stop("the parameters k for shifting the data x series are not set with correct size!")
			#	}
		#}
		#else
		#{
		#	if(length(x)/2/xlen!=length(ks)+1)
		#		{
		#			stop("the parameters k for shifting the data x series are not set with correct size--error02!")
		#		}
		#}
		wmatrix<-mvar
		if(aggregated){
			xRepeat<-1;
		} else{
			xRepeat<-2
			wmatrix<-rep(mvar,2) }
		#ks<-c(1,2,3)
		#ks<-c(0,ks)
		#ks<-c(1,2,3)
		#if(shift.index==1)
		#{
		#	ks<-c(0,ks)
		#} else if(shift.index==length(ks)+1)
		#{
		#	ks<-c(ks,0)
		#}else
		#{
			#insert
		#	ks<-c(ks[1:shift.index-1],0,ks[shift.index:length(ks)])
		#}
		
		ks<-rep(ks,xRepeat)
		#ks<-rep(ks,rep(xlen, length(ks)))
		#cat("before -- xinput:",xinput,"\n")
		xinput<-x+ks
		#cat("after -- xinput:",xinput,"\n")
		if(missing(d))
		{
			d<-pmt_saturated_reading
			if(ylog)
			{
				a<-log(a)
				#d<-log(d)
			}
		}
		if(missing(a))
		{
			a<-pmt_lowest_reading
			if(ylog)
			{
				#a<-log(a)
				d<-log(d)
			}
		}

		#pred<-a+(d-a)/(1+exp((xmid-xinput)/scal))^g
		pred<-f5pl(c(a,d,xmid, scal,g),xinput);
		m<-weight.matrix(pred,model.weight,order)
		if(!missing(mvar))
		{
			#cat("duplicated wmatrix:",wmatrix,"\n")
			m<-wmatrix;
		}
		#cat("\ny:",y, "\n")
		#cat("\npred:", pred, "\n")
		#cat("ks:", ks, "\n")
		#cat("RSS:", sum((log(y)-log(pred))/m*(log(y)-log(pred))/m),"\n")
		(log(y)-log(pred))/m #sqrt(abs(pred))#pred #^2#*(pred)
	}
	
	#this is the function used to calcuate the weight matrix 
	#based on the prediction 
	#model: uniform, 1 
	#		sqrt, sqrt(pred)
	#		log, log(pred)
	#		exp, (pred)^order  order=1,2,3,4 
	weight.matrix<-function(pred, model="uniform", order=1)
	{
		m<-switch(model,
			"uniform"=rep(1,length(pred)),
			"sqrt"=sqrt(abs(pred)),
			"log"=log(abs(pred)),
			"power"=(pred)^order,
			"exp"=sqrt(1/exp(pred)),
			{print("WARNING: unknow model chosen for the weight matrix of 5PL fitting. Using the default one.")
				pred}
		)
		m
	}
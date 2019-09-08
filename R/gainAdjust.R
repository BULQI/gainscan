##=======gainAdjust.R==========
#
## This is the module to do the gainAdjust. The job is to based on the
## readings from multiple gain setting, "summarize" the reading in order
## to 1)minimize the saturation of the high reading points
## as well as 2)increase the low reading points.
##
## The model is based on a 5 parameter logistic function.
## 
## We assume all the readings should be well described a 
## common parameterized logistic function. The different
## protein concentration as well as the different affinity
## together determine/shift the points horizontally
## (on the x-axix). We just need to assume one
## function with count set of parameters (5 of them)
## and but also "align" by fitting to estimate 
## the shift/offset on the x axis
##
## Another thing to note is that we assume the maximum and 
## minimum values of the function is set 0 and 65535 respectively.
## We could allow them to vary and estimate it from the data, but
## proctically this fitting is not ideal. Fixing them makes the 
## fitting more reasonable and easy to converge. (***TODO: allow the   <==============
## maximum one fixed?? allowing minimum value to be varying??Need to
## try!!)
##	a=0.01, d=65535
##
## another thing to mention: the form of the logistic function we assume
## here is the one assume a log transformed x input, since it takes an exponetial
## in the formula. For y, it doesn't have to be log transformed. Either way, it should work,
## but it is said that log transformed is more symetric!!!
## 	We will try both anyway
##
###logistic function 5pl

#for the 
###model 5pl  --- deprecated ---obsolete
##  y<-a+(d-a)/(1+(c/x)^b)^g
### y - log values
### a - smallest y log
### d - hightest y log
### c - might point of x
### b - slope
### g - the asymetric factor
#####

###updated model now, is
###model 5pl
##  y<-a+(d-a)/(1+exp((xmid-x)/scal))^g
### y - log values
### a - smallest y log
### d - hightest y log
### xmid - middle point of x
### scal - slope at the middle point
### g - the asymetric factor
#####
##=======update 
###8/7/2017, move all the residual function to another file residualFunc.R
##
###========update
##1/24/2019, now change the model from 5pl to 
## exponential with background as
## see Tom's note about the
##    log(Y)=log(B+phi*v^d)+eps
##         Y=2^16-1 if(Y>2^16-1)
####======================================================


#------------start the code------------

##first define constants
#Maximizing the Dynamic Range of a Scan
#The complete dynamic range of the scanner is being used when you see a range of intensities on the image from 1 to 65535. 
#A pixel with an intensity of 65535 is "saturated". Saturated pixels represent a condition in which there are more photons
#detected than the PMT can process, or the output of the PMTs exceeds the range of the analog-to-digital converters. 
#A saturated pixel is not an accurate measurement of the signal from the pixel, so it is imperative to 
#set the PMT Gain to avoid saturation.
#from GenePix_6.0_manual.pdf" page 40
#the scanner is using 16 bit encoding, so we have 2^16-1 as the maximum value.
pmt_saturated_reading<-65535;
pmt_lowest_reading<-1
a<-pmt_lowest_reading;
d<-(pmt_saturated_reading-5);

#'@include residualFunc.R


#this function is used to prepare the input of x
# in order to call the fitting function,
# since when we call nlsLM, it check for the input y and x,
# it will complain, when lengths of y and x are not identical
# y, dependent variable, a data matrix or data frame, each row is for one set of data
# x, independent variable, a vector of length of ncol of y

# pos, the position of unchanged/reference xs, 1 by default, 
##		this one is used ONLY when we want to move the reference series to the beginning
##		of the data. NOw this is not necessary. The code in the other
##		part is taking in the ref.index for the reference series.

# return value, a list of two vectors of equal length, y and x

##description, in this function, we prepare the input, transform the
## input into vectors (both y and x), and also based on "pos",
## we moved the referenced y and x into the first slot
## so then in the fitting module, we can assume this and 
## add the k, the shifting parameters, to the following 
## xs. 
### y must be either dataframe or matrix in a format of genes by pmts
### x could be a vector indicating one set of pmts and
###			or it could be matrix or dataframe has same dimension of y
### we also do estimation of variance
#'@export
gainAdjust.prepareInput<-function(y, x, 
		#pos=1, 
		var.log=F
		)
	{
		#need to make sure ncol of y equals to length(x)
		if(class(y)=="numeric") #this is the one set of data, so return everything
		{
			#check for correct input
			if(length(y)!=length(x))
			{
				stop("the input y and x don't have correct format/length, please check!")
			}
			return(list("y"=y,"x"=x, "var"=1))
		}
		
		if(class(y)!="data.frame" && class(y)!="matrix")
		{
			stop("unsupported data format for input y!")
		}
		
	#now we are dealing with data matrix of frame
		
		#check for correct input
		if(class(x)=="matrix"|class(x)=="data.frame")
		{
			x<-as.matrix(x)
			if(dim(y)[2]!=dim(x)[2]|dim(y)[1]!=dim(x)[1])
			{
				stop("ERROR: the input y and x don't have correct format/length, please check!")
			}
		} else if(class(x)=="numeric")
		{
			if(dim(y)[2]!=length(x))
			{
				stop("ERROR: the input y and x don't have correct format/length, please check!")
			}
		}else
		{
			stop("unsupported data format for input x!")
		}
		
		#doing x first
		x_t<-c();
		if(class(x)=="numeric")
		{
			x_t<-rep(x,dim(y)[1])
		} else
		{
			x_t<-rep(0,dim(x)[1]*dim(x)[2])
			for(i in c(1:dim(x)[1]))
			{
				x_t[c(1:dim(x)[2])+(i-1)*dim(x)[2]]<-x[i,]
				
			}
		}
		
		#now doing y
		y_t<-rep(0,dim(y)[1]*dim(y)[2])
		var_t<-rep(0,dim(y)[1]*dim(y)[2]/2)	
		
		for(i in c(1:(dim(y)[1])))
		{
			y_t[c(1:dim(y)[2])+(i-1)*dim(y)[2]]<-y[i,]
			#y_t[c(1:dim(y)[2])+((i)*2-1)*dim(y)[2]]<-y[i*2,]
			#if(var.log)
			#{
			#	var_t[c(1:dim(y)[2])+((i-1))*dim(y)[2]]<-(log(y[i*2-1,]/y[i*2,]))^2
			#}else {
			#	var_t[c(1:dim(y)[2])+((i-1))*dim(y)[2]]<-(y[i*2-1,]-y[i*2,])*(y[i*2-1,]-y[i*2,])
			#}
		}
		
		#check first
		if(floor(dim(y)[1]/2)*2!=dim(y)[1]){#so the number of rows are not even, it might be wrong
			warning("Warning: the number of row of the y input matrix is not even, calling doing aggregration might be wrong!!!")
		}  
		
		for(i in c(1:(dim(y)[1]/2)))
		{
			#y_t[c(1:dim(y)[2])+((i)*2-1-1)*dim(y)[2]]<-y[i*2-1,]
			#y_t[c(1:dim(y)[2])+((i)*2-1)*dim(y)[2]]<-y[i*2,]
			if(var.log)
			{
				var_t[c(1:dim(y)[2])+((i-1))*dim(y)[2]]<-(log(y[i*2-1,]/y[i*2,]))^2
			}else {
				var_t[c(1:dim(y)[2])+((i-1))*dim(y)[2]]<-(y[i*2-1,]-y[i*2,])*(y[i*2-1,]-y[i*2,])
			}
		}
		
		##change the parameters
		#if(pos!=1)
		#{
		#	y1<-y_t[c(1:dim(y)[2])]
		#	y_t[c(1:dim(y)[2])]<-y_t[(pos-1)*dim(y)[2]+c(1:dim(y)[2])]
		#	y_t[(pos-1)*dim(y)[2]+c(1:dim(y)[2])]<-y1
		#	# have to do the same anything on x, if x is not a vector
		#	if(class(x)!="numeric")
		#	{
		#		x1<-x_t[c(1:dim(x)[2])]
		#		x_t[c(1:dim(x)[2])]<-x_t[(pos-1)*dim(x)[2]+c(1:dim(x)[2])]
		#		x_t[(pos-1)*dim(x)[2]+c(1:dim(x)[2])]<-x1
		#	}
		#}
		return (list("y"=y_t, "x"=x_t, "var"=var_t))
	}
	
	#------this is the simpliest way to prepare the inits for shifting
	#---- kind of arbitrary and not working well
	#---Depricated!!!
	#-- It works like this.
	#	--- #arbitrarily choose the first one in the data set as the reference
	#	---           all the lines shifted towards this one
	#	--- #to determine the init k values, we simply compare the biggest Y in 
	#	---				in each data series and take log of the quotient scaled by
	#	---				log(100), another aribitrary choice.
	#The accessary function for preparint intial values for nlsLM
	#this is necessary, because we might have many data series.
	#we need this for actually prepare the initial values
	
	#take in the input, before the prepareInput 
	#and make the initial values in order to do nlsLM fitting
	#it returns an array which is one less the total series (row lengths)
	#Always use the first one (row) as the reference
	#all the following rows are shifted left or right
	#input
	#	y is the input dataframe, has NOT been prepared. it has row as genes and col as pmts
	#	aggregated is indicated whether the data has been aggregated, mean the two repeats has been averaged
	#				by default it is not. So in this case we need to give the two repeat the identical shift, k
	#'@export
	###don't use
	gainAdjust.prepareInits<-function(y, aggregated=TRUE)
	{
		#get dimensions
		rlen<-dim(y)[1]
		clen<-dim(y)[2]
		yTemp<-y#y[seq(1,rlen,by=2),]
		####need to check for the aggregated status, in order to "arrange" data
		if(!aggregated)#need to summarize the data
		{	
			yTemp<-y[seq(1,rlen,by=2),]
			if(rlen%%2!=0)
			{
				stop("the number of array rows are not even, can not do aggregation! Stop!! ");
			}
			for(i in seq(1, rlen, by=2))
			{
				yTemp[(i+1)/2,]=(y[i,]+y[i+1,])
			}
		}
		rlen<-dim(yTemp)[1]
		ref<-yTemp[1,clen]
		
		inits<-rep(0,rlen-1)
		for(i in c(2:rlen))
		{
			inits[i-1]<-log(yTemp[i,clen]/yTemp[1,clen])/log(100)
		}
		inits
	}
	
	#'@title prepare initial values for fitting shifts 
	#'@description generate the initial values for fitting shifts with 
	#'	a model of the 5-parameter logistic function.
	#'@details this is a more complicated way to prepare the initials for shifting. 
	#'	The assumption is that the signal strengths of a feature in protein microarray
	#'	(ProtoArray, invitrogen) acquired with different PMT gains follow 
	#'	a 5-parameter logistic function (5pl) as \cr
	#'	\ifelse{html}{\out{<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-MML-AM_CHTML"></script>
	#'						$$y = a + {{d-a} \over {(1 + e^ {{(xmid - x)} \over scal})}^g} $$
	#'				 	}}{\deqn{y=a+\frac{d-a}{(1+e^(\frac{xmid-x}{scal}))^g}}{ASCII}}
	#'	\cr 
	#'	where a is the upper bound; d is the lower bound); 
	#'	scal is the slope of the linear middle part of the line; 
	#'	xmid is x value at which y equals to (a-d)/2; g is the asymmetric factor.\cr
	#'	On the same array, the 5pl lines for different
	#'	feature proteins should follow a similar pattern, but shifted horizontally left or
	#'	right based on the level of the secondary antibodies bound, which is determined
	#'	by the level of printed feature proteins as well as the interaction between
	#'	the feature protein and the testing antibody/protein. In terms of the
	#'	5pl function parameters, the 5pl lines for different feature proteins are identical
	#'	in a, d, scal and g, but different in xmid; That is to say that if we shift these 
	#'	lines horizontally with the correct amounts, the lines can merge into one line.
	#'	Therefore, we fit these lines together to one 5paremeter logistic function with
	#'	identical a, d, scal and g, as well as one xmid for one reference line and one
	#'	shift for each other lines. In this function, we use a linear model to generate
	#'	the initial values for these shifts parameter of the fitting. \cr
	#'	The algorithm takes a local linear model to generate the initial values like this,
	#'	
	#'	\itemize{
	#'	\item {1. carefully choose the data set as the reference\cr}
	#'		{
	#'		all the lines shifted towards this one.
	#'		the reference line could be anywhere in the data.
	#'		Theoretically, it can
	#'		be any one, but empirically we pick the one line, whose maximum value is closest
	#'		to d/2 (65535/2). This way the reference line has the best coverage over the range of 
	#'		5pl and easy to converge for the fitting.
	#'		}
	#'	\item {2. determine the initial k values\cr}
	#'		{
	#'		we simply compare the biggest Y in 
	#'		in each data series. Two cases are possible \cr
	#'			\itemize{
	#'				\item {1}{ Ymax_i < Ymax_r, nothing to do }
	#'				\item {2}{ Ymax_i >Ymax_r, find the max Yj_i in the series to be smaller then Ymax_r}
	#'			}
	#'		Now we have a pair (Yj_i, Xj_i).
	#'		}
	#'	\item {3.Determine where does this pair belongs to inside the Yr data line\cr}
	#'		{
	#'		Find (Y(k-1)_r, Yk_r), where Y(k-1)_r<Yj_i<Yk_r. Then 
	#'		simply using a linear model to determine Xj_i_pred for Yj_i according
	#'		to reference data . It could possible happen that Yj_i is not inside
	#'		reference range, meaning Yj_i<Ymin_r. In this case, we simply assuming
	#'		Yj_i == Ymin_r, using X(Ymin_r) as the predication and shift the data.
	#'		Hopefully, there will not be many cases like this.
	#'		}
	#'}
	#'to estimate the shift k, we take 
	#'\eqn{k=log(Xj_i/Xj_i_pre)}
	
	#The accessary function for prepare intial values for nlsLM
	#this is necessary, because we might have many data series.
	#we need this for actually prepare the initial values
	
	#take in the input, before the prepareInput 
	#and make the initial values in order to do nlsLM fitting
	#it returns an array which has length identifical to the total series (row lengths)
	#
	#'@param y dataframe this is raw data and has NOT been prepared. 
	#'		it has row as genes and col as pmts
	#'@param x numeric vector It contains the PMT settings for each array. 
	#'@param data.aggregated Parameter indicates whether the data has been 
	#'		aggregated, mean the two repeats has been averaged. 
	#'		By default it is not. In this case, we will generate 
	#'		the identical shift, k, for the two repeats.
	#
	#'@param ref.line list mean to use an external line as the reference line.
	#'		if missing then, use the internal line as the reference.
	#'		It is a list containing two items. The first one is the ref.index 
	#'		for the index of the reference line in the first data set. The
	#'		second is the actually lines. this could be two line for unaggregated
	#'		data or one for aggregated data.
	#'
	#'@return a data list contain two things. First is the ref.index, the index
	#'		of the line in the data set. If we are using the external reference 
	#'		line, then this is -1. The second is estimated initial values for
	#'		shifts. The shifts INCLUDE the shift for the reference line with
	#'		a value of zero. 
	#'		Note: 1) the prepared inits vector always is in a format of aggregated data
	#'		; 2) in log'ed fold; and
	#'		the output doesn't include ref.line when we reference to the external lines
	#'@export
	gainAdjust.prepareInitsLM<-function(y, x, aggregated=TRUE, ref.line)
	{
		#get dimensions
		rlen<-dim(y)[1]
		clen<-dim(y)[2]
		yTemp<-y#y[seq(1,rlen,by=2),]
		####need to check for the aggregated status, in order to "arrange" data
		if(!aggregated)#need to summarize the data
		{	
			yTemp<-y[seq(1,rlen,by=2),]
			if(rlen%%2!=0)
			{
				stop("the number of array rows are not even, can not do aggregation! Stop!! ");
			}
			#for(i in seq(1, rlen, by=2))
			#{
			#	yTemp[(i+1)/2,]=(y[i,]+y[i+1,])/2
			#}
			in.y<-seq(1,rlen,by=2)
			yTemp[(in.y+1)/2,]=(y[in.y,]+y[in.y+1,])/2
			if(!missing(ref.line))
			{
				ref.line$y[1,]=(ref.line$y[1,]+ref.line$y[2,])/2
			}
		}
		rlen<-dim(yTemp)[1]
		#ref<-yTemp[1,clen]
		
		#sort x in case X is not in order
		idx.srtX<-order(x)  #increasing order
		
		inits<-rep(0,rlen)
		Y<-yTemp[,idx.srtX]
		X<-x[idx.srtX]
		idx.ref<-which.min(abs(pmt_saturated_reading/2-Y[,clen]))##here yTemp is the sum of the duplicated data, so we choose the reference to be the one closest to MaxF/2
		#idxs.arr<-yTemp[-1])
		ref<-yTemp[idx.ref,idx.srtX]
		if(!missing(ref.line))
		{
			idx.ref<- -1
			ref<-ref.line$y[1,idx.srtX]
		}
		
		for(i in c(1:rlen))
		{
			if(i==idx.ref) #this is the reference array, skip to the next,
			{
				next
			}
			#do the job to calculate the shifting ks
				#first compare the Ymax with Ymax_r
			#want to find a Yi that is smaller than ref_max
			Yi<-Y[i,]
			
			j<-1;
			for(j in length(Yi):1)
			{
				if(Yi[j]<=ref[clen])
				{
					break; ##done,
				}
				
			}
			Yi_j<-Yi[j]
			Xi_j<-X[j]
			if(Yi_j>ref[clen]) #even the smaller Y in the search one is bigger than the ref max, what to do?? 
			{
				Xi_r_pred<-X[clen]
			}
			else
			{
				#now that, we find a pair (Yi, Xi) to compare with ref data
				#next need to determine, which region it locates in ref data
				k<-1
				for(k in (clen-1):1)
				{
					if(Yi_j>=ref[k])
					{
						break; #we found a good one 
					}
				}
				
				if(Yi_j>= ref[1]) #do linear model
				{
					a<-(ref[k]-ref[k+1])/(X[k]-X[k+1])
					b<-((ref[k]/X[k])-(ref[k+1]/X[k+1]))/(1/X[k]-1/X[k+1])
					Xi_r_pred<-(Yi_j-b)/a
				}
				else  #Yi is smaller than smallest ref
				{
					Xi_r_pred<-X[k] #in this case, k=1
				}
				
			}
			inits[i]<-log(Xi_r_pred/Xi_j)
		}
		list("inits"=inits, "ref.index"=idx.ref)
	}
	
		
	#need to plot as a function
	#y is the original dataframe and has not been "changed"
	#x is the pmt single array
	#model, in this case, we are using 5PL, but with fixed highest end (fix.high), fixed both end (fix.both) and fixed none (fix.none)
	#ylog, in this case, x input is alway natural scale and will always be log transform in here. 
	#But , ylog determines whether y data input has been log transformed. If not transformed (ylog=F), then a, d 
	# inputs will be using natural scaled. The expected y data as determined by the equation will be log transformed
	# in order to do plotting. If data input are log-transformed (ylog=T), then a and d inputs will have to be transformed
	# . But since the data are transformed, meaning the fitted equation is also in log-transformed format, so we don't
	# have to transformed expected y (calculated by equation).
	#'@export
	gainAdjust.plot<-function(LMfit,y, x, ylog=TRUE,aggregated=FALSE,
			model="fix.both",
			a,d, xmid, scal, g, ref.index=1,
			filename=NULL)
	{
		if(missing(LMfit))
		{
			stop("no LM nls fitting object has been specified!")
		}
	
		if(class(LMfit)!="nls.lm"&&class(LMfit)!="list")
		{
			stop("the input is not valid. A nls.lm or list object with par field is needed!")
		}
		#prepare the parameters
		#determine the parameters first
		#define the default case of 3PL
		#xmid<-LMfit$par[1]
		#scal<-LMfit$par[2]
		#g<-LMfit$par[3]
		k_index<-4
		if(missing(a))
			{
				a<-pmt_lowest_reading
			}
		if(missing(d))
			{
				d<-pmt_saturated_reading
			}
		if(ylog)
		{
			a<-log(a)
			d<-log(d)
		}
		
		switch(model,
				"fix.both"={	
						#this has been taken care of as the default case in above
						xmid<-LMfit$par[1]
						scal<-LMfit$par[2]
						g<-LMfit$par[3]
						#do it again here, because we can not leave it empty
						k_index<-4
					},
				"fix.high"={
						#be careful here, 4PL is not the real 4PL, it is 5PL but fixed d, maximum level
						a<-LMfit$par[1]
						xmid<-LMfit$par[2]
						scal<-LMfit$par[3]
						g<-LMfit$par[4]
						k_index<-5
					},
				"fix.none"={
						#all five parameter
						a<-LMfit$par[1]
						d<-LMfit$par[2]
						xmid<-LMfit$par[3]
						scal<-LMfit$par[4]
						g<-LMfit$par[5]
						k_index<-6
				},
				"fix.all"={
					if(missing(xmid)|missing(scal)|missing(g))
					{
						stop("please specify parameters (xmid/scal/g)!!!")
					}
					k_index<-1
					
				},
				stop("unkown model specified")
			)
		#now doing the plot
		if(!is.null(filename))
		{
			jpeg(filename=filename)
		}
		#op<-par(mfrow=c(2,1), mar = c(3,5,2,1))
		#need to figure out the range of y and x
		#for y, we know it is (0.1, 65535)
		
		#for x, we need to add ks
		xmin<-log(min(x))+min(LMfit$par[c(k_index:length(LMfit$par))])
		xmax<-log(max(x))+max(LMfit$par[c(k_index:length(LMfit$par))])
		
		ymin<-min(log(a), min(y))
		if(ylog)
		{
			ymin<-min(a,min(y))
			ymax<-d
		} 
		else ##not log, then we need to turn everything into log, that is how we plot
		{
			ymin<-min(log(a),log(min(y)))
			ymax<-log(d)
		}
		
		plot(c(xmin-0.2,xmax+0.3), c(ymin-0.2,ymax+0.2), bg = "black", cex = 0.5, main="data", type="n")
		#plot(c(xmin,xmax+0.3), c((0.1),(66535)+5000), bg = "black", cex = 0.5, main="data", type="n")
		#points(log(x), log(y[1,]), col="grey", lty=2) 
		params<-LMfit$par[c(k_index:length(LMfit$par))]
		if(ref.index==1){
			params<-c(0,params)
		}else if (ref.index==length(params+1)){
			params<-c(params,0)
		}else{
			params<-c(params[1:(ref.index-1)],0,params[ref.index:length(params)])
		}
		yParams<-y #[-1,]
		if(!aggregated)
		{
			#cat("calling it now for aggregation")
			#yParams<-yParams[-1,];
			#params<-rep(params, rep(2,length(params)));
			yParams<-sqrt(y[seq(1, dim(y)[1],by=2),]*y[seq(2, dim(y)[1],by=2),])#/2
			if(length(params)!=dim(yParams)[1])
			{
				stop("Error: the length of params is not identical to number of y data points!")
			}
			#for(j in c(1:length(params)))
			#{
			#	yParam[j,]<-(y[j*2-1,]+y[j*2,])/2
			#
			#}
			
		}
		for(i in c(1:length(params)))
		{
			#if(i!=304)
			#{
			# next;
			#}
			#cat("\ny:", yParams[i,], "\n")
			#cat("param[i]:", params[i], "\n")
			#points(log(x)+params[i], log(yParams[i,]), col="grey", lty=2)
				if(ylog) #this means that the data is logged, so don't do log again
				{
					lines(log(x)+params[i], yParams[i,], col=i, lty=2) 
				}
				else
				{
					#------>
					#cat("-->long:",log(yParams[i,]), "\n" )
					lines(log(x)+params[i], log(yParams[i,]), col=i, lty=2) 
				}
		}
			
		#now plot the fitted value
							
		xx<-seq(xmin-0.1,xmax+0.2, by=(xmax+0.2-xmin+0.1)/1000)
		#cat(xx,"\n")
		yy<-f5pl(c(a,d,xmid,scal,g),xx);#yy<-a+(d-a)/(1+exp((xmid-xx)/scal))^(g)  #<--changed now
		#cat(yy, "\n")
		cat("a:", a, "\td:", d, "\txmid:", xmid, "\tscal:",scal, "\tg:", g, "\n")
		if(!ylog)
		{
			#---->
			yy<-log(yy)
		}
		lines(xx,yy,col=2, lty=1,lwd=2)
		#--plot(1:(LMfit$niter+1), log(LMfit$rsstrace), type="b",
		#--main="log residual sum of squares vs. iteration number",
		#--xlab="iteration", ylab="log residual sum of squares", pch=21,bg=2) 
		#par(op)
		if(!is.null(filename))
		{
			dev.off();
		}
	}
	
	
	#'@title aggregate the duplicated data, either by arithmatic or geometric mean
	#'@description we assume the input data eitehr data matrix or an EListRaw data \code{\link{EListRaw-class}}.
	#'	In either case, we assume the data are duplicated and arranged in the adjacent lines.
	#'@details We can do either "geometric" or "arithmetic" mean between the duplicated data. 
	#'	Also for the EListRaw data object, we will aggregate the E, C, cgenes and genes data tables.
	#'	but we don't change the background data (neither Eb nor Cb). 
	#'@param y data object, either data matrix or EListRaw.
	#'@param mode string to indicate whether to do geometric or arithmetic mean.
	#'@return aggregated data
	#'@seealso \code{\link{EListRaw-class}}
	#'@export
	gainAdjust.aggregate<-function(y,  mode="geometric")
	{
		if(missing(y))
		{
			stop("data missing, please specify")
		}
		y.aggregated<-c(); #set up output
		
		switch(class(y),
		"matrix"={
			nrow<-dim(y)[1]
			#ncol<-dim(y)[2]
			
			#assuming the duplicated data are next to each other
			if(nrow!=floor(nrow/2)*2)
			{
				stop("data is corrupted. the number of data rows is not even!")
			}
			idx<-seq(1, nrow, by=2)
			y.aggregated<-y[idx,]
			if(mode=="geometric") #this means that the data is logged, so don't do log again
				{
					y.aggregated<-sqrt(y.aggregated)*sqrt(y[idx+1,]) 
				}
			else if(mode=="arithmetic")
				{
					#------>
					#cat("-->long:",log(yParams[i,]), "\n" )
					y.aggregated<-(y.aggregated+y[idx+1,])/2
				}
			else
				{
					stop("unknow mode, specify either arithmatic or geometric")
				}
		
		}, 
		"EListRaw"={ ##need to take care both the datamatrix and gene dataframe
			y.aggregated<-y;

			#do E data first
			nrow<-dim(y$E)[1]
			#ncol<-dim(y)[2]
			
			#assuming the duplicated data are next to each other
			if(nrow!=floor(nrow/2)*2)
			{
				stop("data is corrupted. the number of data rows is not even!")
			}
			idx<-seq(1, nrow, by=2)
			y.aggregated$E<-y$E[idx,]
			#storage.mode(y.aggregated$E)<-"double"
			y.aggregated$genes<-y$genes[idx,];
			if(mode=="geometric") #this means that the data is logged, so don't do log again
				{
					y.aggregated$E<-sqrt(y.aggregated$E)*sqrt(y$E[idx+1,]) 
				}
			else if(mode=="arithmatic")
				{
					#------>
					#cat("-->long:",log(yParams[i,]), "\n" )
					y.aggregated$E<-(y.aggregated$E+y$E[idx+1,])/2
				}
			else
				{
					stop("unknow mode, specify either arithmatic or geometric")
				}
			
			#doing the control proteins
			nrow<-dim(y$C)[1];
			#assuming the duplicated data are next to each other
			if(nrow!=floor(nrow/2)*2)
			{
				stop("data is corrupted. the number of data rows is not even!")
			}
			idx<-seq(1, nrow, by=2)
			y.aggregated$C<-y$C[idx,]
			y.aggregated$cgenes<-y$cgenes[idx,];
			if(mode=="geometric") #this means that the data is logged, so don't do log again
				{
					y.aggregated$C<-sqrt(y.aggregated$C)*sqrt(y$C[idx+1,]); 
				}
			else if(mode=="arithmatic")
				{
					#------>
					#cat("-->long:",log(yParams[i,]), "\n" )
					y.aggregated$C<-(y.aggregated$C+y$C[idx+1,])/2
				}
			else
				{
					stop("unknow mode, specify either arithmatic or geometric")
				}
			#done, we don't other 
		}, #case two
		{
			stop("Not correct data format. data matrix is expected!")
		} #default case
		)
		
		y.aggregated
	}
	
	#this is the function to align the data series into one
	#5pl lines based on the fiting with shifted parameters
	#the input is 
	# LMfit, with the shifted parameters
	#	y in datafame format
	#	x the original x data series in PMTs
	#return value
	#	dataframe for shifted xs, in order of y . And 
	#   this is not log transformed
	#'@export
	gainAdjust.alignData<-function(LMfit,x, aggregated=F,model="fix.both", ref.index=1)
	{
		k_index<-4
		
		switch(model,
				"fix.both"={	
						#do it again here, because we can not leave it empty
						k_index<-4
					},
				"fix.high"={
						#be careful here, 4PL is not the real 4PL, it is 5PL but fixed d, maximum level
						k_index<-5
					},
				"fix.none"={
						#all five parameter
						k_index<-6
				},
				"fix.all"={
						k_index<-1
				},
				stop("unkown model specified")
			)
		params<-c(LMfit$par[c(k_index:length(LMfit$par))])
		if(ref.index==1){
			params<-c(0,params)
		}else if(ref.index==length(params)+1){
			params<-c(params,0)
		}else{
			params<-c(params[1:ref.index-1], 0, params[ref.index:length(params)])
		}
		
		
		if(!aggregated)
		{	
			params<-rep(params, rep(2,length(params)));
				
		}
		params<-rep(params,rep(length(x),length(params)))
		
		params<-matrix(params,ncol=length(x),byrow=T)
		xFrame<-data.frame(matrix(rep(x,dim(params)[1]),ncol=length(x),byrow=T))
		xFrame*exp(params)
		
	}
	#'@title fit the shift for single feature with 5-parameter logistic function
	#'@description fit the logarithmic line shift for individual feature based 
	#'	on predetermined 5-parameter logistic function (5pl).
	#'@details In this function the single feature data are fitted with 5pl to
	#'	estimate the shift, k. Here what it means is that the single feature data
	#'	should follow a 5pl with fixed parameters of a, d, scal and g as to a 
	#'	reference line. The feature line and the reference will merge into one
	#'	line if we keep the reference line fixed and move/shift the feature line
	#'	horizontally by k. 
	#'@seealso \code{\link{gainAdjust.prepareInitsLM}} 
	#'	\code{\link{gainAdjust.fit5plArray}}
	#'@param ydata data matrix holding the raw data to be fitted with 
	#'	features in row and PMT gain settings in column
	#'@param x numeric vector the PMT gain settings to acquire the data.
	#'	it has the order as to that of column of ydata
	#'@param par.5pl nlm fitting object holding the parameters of 5pl
	#'@param data.aggregated logic indicating whether the data have been aggregated
	#'@param debug logic when set to be true, the debugging plot will be shown
	#'@param ref.index numeric an index used to indicate whether the reference line
	#'	is inside the initsLM (see next input parameter)
	#'@param initsLM a numeric vector contains the initial values for the line shift k. please
	#'	\code{\link{gainAdjust.prepareInitsLM}}. Want to emphasize that both parameter ref.index
	#'	and initsLM are the output of the function \code{\link{gainAdjust.prepareInitsLM}}. 
	#'	initsLM contains the shift for the reference line with a value of zero.
	#'
	#'@return a data list contains one field "par". This is the par estimated for the
	#'	shift k for each feature line. For the unaggregated data, the shift k parameters
	#'	for duplicated data are identical. Also, for a internal reference line (ref.index!=-1)
	#'	, the output doesn't contain shift for it. Otherwise (ref.index=-1), it contain all
	#'	shifts for the data feature lines. The shifts are in a logged fold format.
	#'@export
	fitShifts5plSingle<-function(ydata, x, par.5pl, data.aggregated=F, debug=F, ref.index, initsLM)
	{
		if(missing(ydata))
		{
			stop("please specify the expression data!")
		}
		if(missing(x))
		{
			stop("please specify the PMT setting data!")
		}
		if(missing(par.5pl))
		{
			stop("please specify the fitting parameters!")
		}
		#if(length(var)!=1&&length(var)!=dim(ydata)[1])
		#{
		#	stop("the var is not set up correctly")
		#}
		
		if(missing(ref.index))
		{
			stop("***no input for ref.index, so run to find it again...\n")
			#initsLM<-gainAdjust.prepareInitsLM(ydata,x, aggregated=data.aggregated)
			#ref.index<-initsLM$ref.index
		}
		if(missing(initsLM))
		{
			stop("***no input for initial values, please specify")
		}
		if(class(initsLM)!="numeric")
		{
			stop("**unknown format for initLM, please check.")
		}
			#parStart<-c(6,0.1,0.13,-0.2,-0.5,-0.55, -0.5,1)
				
			#-gainAdjust.fitFinal<-nls.lm(par=initsLM$inits[-initsLM$ref.index], 
			#-	fn=gainAdjust.fnResShift,
			#-	#y=log(yinput.aligned), 
			#-	y=(yinput), 
			#-	x=log(xinput), xlen=xlen,aggregated=data.aggregated,
			#-	ylog=data.ylog, model.weight="log", #in here, this fnResShift doesn't take any model, it hardcodes using log transform
			#-	order=1,d=(d),a=(a), 
			#-	xmid=gainAdjust.fit$par[1],scal=gainAdjust.fit$par[2],g=gainAdjust.fit$par[3],
			#-	shift.index=initsLM$ref.index,
			#-	control = nls.lm.control(nprint=1, maxit=100)
			#-	)
		#the statement block for running fitting with inividual data sets	
		#{
		
			fitShift<-list();
			fitShift$par<-initsLM;#[-initsLM$ref.index]; #initialize the output
			#--->xsim<-seq(4,7,by=0.01)
			#--->plot(exp(xsim), f5pl(c(a,d, gainAdjust.fit$par[1:3]),xsim),col=2, lwd=2, type="l", log="xy", lty=2)
			cat("Fitting the shifting parameters......\n")
			nprint<-0;
			if(debug)
			{
				nprint<-1;
			}
			cat("\nStart fitting...");
			cnt<-0; #indicator for printing progresses
			
			for(i in 1:length(initsLM))
			{
				if(debug){
					cat(i,"/", length(initsLM)
					,"...\n")
				}
				if(i==ref.index)
				{
					fitShift$par[i]<-0
					next;
				}
				xlen<-length(x)
				#yvar<-(log(ydata[i*2-1,]/ydata[i*2,]))^2
				#idx.zero<-which(yvar==0)
				#if(length(idx.zero)>0)
				#{
				#	yvar[idx.zero]<-min(yvar[-idx.zero]) ##to get rid of zero variance.
				#}
				#yvar<-1
				#do fitting individually
				####now based on data aggregation we need to do fitting
				gainAdjust.fitS<-NULL
				if(!data.aggregated){
					
				gainAdjust.fitS<-nls.lm(par=fitShift$par[i], fn=gainAdjust.fnResShiftSingle,
					#y=log(yinput.aligned), 
					y=c(ydata[i*2-1,], ydata[i*2,]), 
					x=log(x), #xlen=9,
					aggregated=data.aggregated,
					#ylog=ylog, model.weight="power", #in here, this fnResShift doesn't take any model, it hardcodes using log transform
					#order=1,
					a=par.5pl[1],d=par.5pl[2], 
					xmid=par.5pl[3],scal=par.5pl[4],g=par.5pl[5], mvar=1,#yvar,
					#shift.index=initsLM$ref.index,
					control = nls.lm.control(nprint=nprint, maxit=100))
				
				} else {
					gainAdjust.fitS<-nls.lm(par=fitShift$par[i], fn=gainAdjust.fnResShiftSingle,
					#y=log(yinput.aligned), 
					y=c(ydata[i]), 
					x=log(x), #xlen=9,
					aggregated=data.aggregated,
					#ylog=ylog, model.weight="power", #in here, this fnResShift doesn't take any model, it hardcodes using log transform
					#order=1,
					a=par.5pl[1],d=par.5pl[2], 
					xmid=par.5pl[3],scal=par.5pl[4],g=par.5pl[5], mvar=1,#yvar,
					#shift.index=initsLM$ref.index,
					control = nls.lm.control(nprint=nprint, maxit=100))
				}
				if(!debug){
					
					if(ceiling((i*100)/length(initsLM))>=cnt*5)
					{
						#cat("fitting shift:")
						cat(paste(cnt*5,"%...",sep=""))
						#cat(gainAdjust.fitS$par[1])
						#cat("\n")
						cnt<-cnt+1;
					}
					
				}	
					 #
				if(debug){
					if(!data.aggregated){
					lines(exp(log(x)+gainAdjust.fitS$par[1]),ydata[i*2,],col=i)
					lines(exp(log(x)+gainAdjust.fitS$par[1]),ydata[i*2-1,],col=i)
					} else {
					lines(exp(log(x)+gainAdjust.fitS$par[1]),ydata[i,],col=i)
					}
				}
				
				fitShift$par[i]<-gainAdjust.fitS$par[1];
				flush.console();
			}#end of for loop
		#}
		cat("Done!!!!\n")
		if(ref.index!=-1)
		{
			fitShift$par<-fitShift$par[-ref.index]
		}
		fitShift
	}
	
	#'@title fit 5-parameter logistic function with one array data
	#'@description fit the 5-parameter logistic function with one single array data. 
	#'	Both feature data and control data will be fitted and adjusted for the input
	#'	single array data.
	#'@details 
	###------------>updated 8/4/2017, feng
	#' the function to fit 5PL for parameters and shifts for INDIVIDUAL array with multiple PMT readings
	#' In this function, we first 1)do fitting with fixed a and d to estimate xmid scal, g as well 
	#'	all shifts with uniform weights.
	#'then we 2) do fitting with the fixed 5PL parameters to estimate individual shifts of x's under log-ed data.
	#'Then in the third step, 3) estimate 5PL parameters, a, xmid, scal and g (fixing d) with shifted data
	#' under log weighted matrix.\cr
	#'	If in the last step, the fitted a is negative, we just simply go back to the original fitted 5PL. 
	### NOTE: somehow, the first fitting fixing both has to be uniform weights to get the best fitting.
	#'	also want to mention, the threshold is important, since it is hard to converge for fitting
	#'	low signal feature data. Now incremental thresholds are implemented to make it work by removing
	#'	problematic ones.
	##8******TODO, get a better init for parameters in the 1st 3p fitting --DONE!!!
	###Input: 
	#'@param y data list is the list containing expression and control expression
	#'\itemize{
	#'		\item{y$C the data matrix for control expression in array with different PMT readings} 
	#'		\item{y$E the data matrix for row is for target expression in array with different PMT readings}
	#'		\item{y$gene  meta data for gene} 
	#'		\item{y$cgene meta data for control gene}
	#'}
	#'@param x the vector containing the PMT gain. In the fitting, we always want to do log(xdata)
	#'			elements in x must be arranged in an order identical to the column order 
	#'			of PMT settings for array in
	#'			
	#'@param weight.mode char string for advanced users. modify them if you know what you are doing.
	#'@param data.ylog logic indicating whether y data are log-transformed
	#'@param data.aggregated logic indicating whether y data are aggregated, 
	#'	meaning the two replicates are averaged arithmatic and geometric.
	#'			(the mode of aggregation was done by user)
	#'@param fit.aggregated and aggregated.mode are describing how to do fitting 
	#'			either using aggregated or replicated data.  Default is raw non-aggregated data. 
	#'			The mode is either "arithmetic" and "geometric". 
	#'@param fit.mode string to specify whether to use control ("control") or 
	#'	control+Expression data ("all"). By default, only
	#'			fit the control data. The concern is the running time. 
	#'	Obviously including more data takes too long to finish.
	#'@param block.size numeric the number of blocks to run fitting. The total numer of blocks is 48. 
	#'			Ideally we should run fitting with one block, but sometimes the features in one block 
	#'	are not enough to fit a good 5pl. The reason biggest size is 12 control blocks, 
	#'	with about 1000 points (500 unique points), which takes not too long to finish either.
	#'@param threshold, the lowest value that maximum PMT gain setting expression should achieve. If data
	#'			points with all expression level below this threshold are not included for fitting 5PL
	#
	#'@param arrayName the name of the array. this is not necessary, by help to write debugging plots names
	#'@param order, is only used by the residual function of the nlm.fit. is used for generation of exponential weight matrix, the order of the power. Not used if other type of weight matrix is applied.
	#'@param PMT.gain data list contain the best PMT.gain for each feature determined by
	#'	5pl fitting. it has two fields, E and C, for feature data and control data respectively.
	#'	when this is null, means we need to determine them in this fitting. Otherwise, we use 
	#'	the input for data adjust.
	#'@param a and d the two parameters for 5pl fitting. 
	#'@return a data list contains following fields, E, E.adj, C, c.adj and PMT.gain.
	#'@seealso \code{\link{gainAdjust.prepareInitsLM}} \code{\link{fitShifts5plSingle}}
	#'	\code{\link{adjust.Matrix}}
	#'@export
	gainAdjust.fit5plArray<-function(y, x, data.ylog=F, data.aggregated=F,
		fit.mode="control",
		fit.aggregated=F, aggregated.mode="geometric", 
		weight.mode, order=1,
		block.size=12,a=1.1, d=65535,threshold=100,
		PMT.gain=NULL,
		arrayName=NULL,
		debug=F
		)
	{
		if(missing(x))
		{
			stop("input \"x\" missing, please specify!");
		} 
		if(missing(y))
		{
			stop("input \"y\" missing, please specify!");
		}
		if(class(y)!="list"){
			stop("ERROR: the input data of y is not in a correct format. List with expression data needed!!")
		}
		if(missing(weight.mode))
		{
			weight.mode<-c("uniform", "log", "sqrt");
		}
		#check for the correct weight mode
		for(i in 1:length(weight.mode))
		{
			if(!is.element(weight.mode[i], c("log", "uniform", "power","exp", "sqrt")))
			{
				stop("Unknow weight mode specified. Only log/uniform/power/exp supported!");
			}
		}
		#need to get data in below based on input,
		blkNum<-max(y$gene[,"Block"])
		if(block.size<1||block.size>blkNum)
		{
			stop(paste("only a block size between 1 and ", blkNum,sep=""))
		}
		if(is.null(arrayName))
		{
			arrayName<-as.numeric(format(Sys.time(), "%OS3"))*1000 ###in milliseconds
			arrayName<-paste("a",arrayName, sep="")
		}
		#prepare the output
		fr<-list();#this is the fitting results
		#dr<-list("E"=y$E[,1], "C"=y$C[,1], "E.adj"=y$E[,1], "C.adj"=y$C[,1], 
		#		"E.gain"=y$E[,1], "C.gain"=y$C[,1]); #this is the gain adjust data output, E/C
															#is the array holding the optimal data for this 
															#array, and PMT.gain is the gain setting
															#for the data
		fit<-vector(mode="list", length=ceiling(blkNum/block.size))
		cshifts<-vector(mode="list", length=ceiling(blkNum/block.size))
		shifts<-vector(mode="list", length=ceiling(blkNum/block.size))
		###now good ?? start doing the fitting, based on the control and block size
		
		ref.line<-list("y"=matrix(0,nrow=2,ncol=length(x),byrow=T),"ref.index"=-1) ##<-this one is used to carry the reference line
		for(j in 1:ceiling(blkNum/block.size)){
			cat("Fitting data by blocks: ",j, "/",blkNum/block.size, "...\n")
			#for the array, we need to get the protein data by block, 
			bg<-(j-1)*block.size+1
			ed<-(j-1)*block.size+block.size
			#cat(paste("bg:",bg, "; ed:", ed,"\n", sep=""));
			if(ed>blkNum)
			{
				ed<-blkNum
			}
			
			idx.prC<-which(is.element(y$cgene[,"Block"],c(bg:ed)))
			ydata<-y$C[idx.prC,]
			if(fit.mode=="all")
			{
				idx.prE<-which(is.element(y$gene[,"Block"],c(bg:ed)))
				ydata<-rbind(y$E[idx.prE,], ydata)
			}
			
			#for the next section, we need to do while loop to take care of the cases where the fitting is not 
			#working properly for some small values series. In these cases, we will keep increasing 
			#threshold to get rid of these data series from the initial fitting.
			flag_threshold_reset<-T
			current_threshold<-threshold
			initsLM<-NULL; #get this ready for using after the while loops
			while(flag_threshold_reset)
			{
			#cat(paste("threshold:",threshold,"\n",sep="") );
			#cat(paste("current_threshold:",current_threshold,"\n",sep="") );
			#now data is ready, rm the low expressed protein, since they are confusing the fitting
			ydata5PL_List<-rmLows(ydata, threshold=current_threshold, index=which.max(x), aggregated=data.aggregated)
			ydata5PL<-ydata5PL_List$m
			xlen<-length(x)
			
			#now we need to do aggregation
			if(!data.aggregated)
			{
				if(fit.aggregated)
				{
					#do aggregation on data
					if(data.ylog)
					{	
						if(aggregated.mode=="geometric")
						{
							cat("WARNING: call to run aggregation on log transformed data. Could lead to NaN values")
						}
					}
					data.aggregated<-T
					ydata5PL<-gainAdjust.aggregate(ydata5PL, aggregated.mode)
				}
			}else  ##input data already aggregated, then do nothing
			{
				#in this case, we ignore whatever "fit.aggregated" says, since there are no way to revert data back
				#so set data.aggregated as true
				data.aggregated<-T; #redundant
			}
			
			#get the initial values for fitting, by picking the data series within the middle range 
			if(j==1)##this is the first block, we need to set up reference line
			{
				initsLM<-gainAdjust.prepareInitsLM(ydata5PL,x, aggregated=data.aggregated)
				if(!data.aggregated){
					ref.line$y<-ydata5PL[c(initsLM$ref.index*2-1,initsLM$ref.index*2),]
				} else {
					ref.line$y<-ydata5PL[initsLM$ref.index,]
				}
				ref.line$ref.index<-initsLM$ref.index
			} else { #in this case, we need to use the reference from block 1 and 
						#also insert this one to do fitting, at location of index 1
				#keep the ref.line unchanged
				initsLM<-gainAdjust.prepareInitsLM(ydata5PL, x, aggregated=data.aggregated, ref.line=ref.line)
				if(!data.aggregated){
					ydata5PL<-rbind(ref.line$y, ydata5PL)
				} else {
					ydata5PL<-rbind(ref.line$y[1,], ydata5PL)
				}
				
				initsLM$ref.index<-1
				initsLM$inits<-c(0,initsLM$inits)
			}
			#parStart<-c(6,0.1,0.13,-0.2,-0.5,-0.55, -0.5,1)
				
			#for three parameters, xmid, scal and g
			parStart<-c(6.5,0.005,0.015,initsLM$inits[-initsLM$ref.index]) #<---3pl, keeping a and d constant
			#data.ylog<-F
			if(data.ylog)
			{
				a<-log(a)
				d<-log(d)
				parStart<-c(6.3,0.1,0.1,initsLM$inits[-initsLM$ref.index]) 
			}
			
			#now we have the data, reay, aggregated
			input<-gainAdjust.prepareInput(ydata5PL,x)
			yinput<-input$y
			xinput<-input$x
			
			cat("\tfitting for shifts and parameters.....\n")
			cat("\t\tcurrent_threshold value: ");
			cat(current_threshold);
			cat(".....\n")
			gainAdjust.fit<-nls.lm(par=parStart, 
				#fn=gainAdjust.fnRes4pShift,  #<----4pl, varying a
				fn=gainAdjust.fnRes3pShift, #<----3pL, keeping a and d constant
				#y=log(yinput), 
				y=yinput,
				x=log(xinput), xlen=xlen, aggregated=data.aggregated, ylog=data.ylog,a=a,d=d,
				shift.index=initsLM$ref.index,
				model.weight="uniform", order=order, control = nls.lm.control(nprint=1,maxiter=100)
				)
			#cat("----doing before calling try catch.........\n")
			flag_threshold_reset<-
			tryCatch(
				{
					summary(gainAdjust.fit)
					#cat("here in tryCatch first block")
					FALSE
				}, #this will break for error if something wrong for the fitting,
				error=function(eMsg){
					#reset the flag and change the threshold
					#flag_threshold_reset<-T
					print(eMsg)
					threshold<-threshold+10
					cat("****rerun the fitting with a new threshold!\n")
					
					TRUE
				},
				warning=function(wMsg){
					#it might never happen???
					print(wMsg);
					cat("\n")
					#cat("executed too")
					#go ahead for now
					#flag_threshold_reset<-F
					FALSE #in case of warning, just go ahead
				},
				finally={
					#do nothing.
					#cat("in here for finllay.")
				}
			
			)#end of try catch block
			if(flag_threshold_reset)
			{
				current_threshold<-current_threshold+10;
			}
			
			#cat("the flag_threshold_reset: ")
			#cat(flag_threshold_reset)
			#cat("\n")
			#cat("threshold value: ")
			#cat(threshold)
			#cat("current_threshold value: ")
			#cat(current_threshold)
			#cat("\n")
			}#----->end of while loop for resetting the threshold.
			if(debug){	
				#sink("debug.txt", append=T)
				#cat("summary of block", j, ":\n")
				#summary(gainAdjust.fit)
				#sink();
			
			#gainAdjust.plot(gainAdjust.fit,log(y), x,ylog , ref.index=initsLM$ref.index)
			#gainAdjust.plot(gainAdjust.fit,y, x,ylog)
			gainAdjust.plot(gainAdjust.fit,ydata5PL, x, ylog=data.ylog , 
								ref.index=initsLM$ref.index,a=a, d=d,aggregated=data.aggregated,
								model="fix.both"
								#,filename=paste("fit1_blk_",j,".jpg",sep="")
							)
			#plot 
			x.aligned<-gainAdjust.alignData(gainAdjust.fit,x, model="fix.both", ref.index=initsLM$ref.index)##<---for 4Pl, 
			input.aligned<-gainAdjust.prepareInput(ydata5PL,x.aligned)
			yinput.aligned<-input.aligned$y
			xinput.aligned<-input.aligned$x
			jpeg(filename=paste(arrayName,"fit1_dot_blk_",j, "%d.jpeg",sep=""))
			plot((xinput.aligned), (yinput.aligned), type="p", main="5-p logistic fitting of intensity vs. PMT gain",
				xlab="PMT gain (log)",ylab="intensity (log)",log="xy"
				)
			#plot(log(xinput.aligned), (yinput.aligned), type="p")
			xsim<-seq(4,7,by=0.01)
			
			lines(exp(xsim), (f5pl(c(a,d, gainAdjust.fit$par[1:3]),xsim)),lty=1,col=3,lwd=2)
			dev.off()
			}
		#done for first run
			
		#now do second fit, only shift data
			#do aggregation now
			#ydata_ag<-gainAdjust.aggregate(ydata, mode="geometric")
			#aggregated_3rd<-T
			#--input<-gainAdjust.prepareInput(ydata,x) #now including even the low expressing data points
			#--yinput<-input$y
			#--xinput<-input$x
			cat("\trefining fitting for shifts .....\n")	
		#get the initial values for fitting, by picking the data series within the middle range 
			initsLM<-gainAdjust.prepareInitsLM(ydata,x, aggregated=data.aggregated)
			#parStart<-c(6,0.1,0.13,-0.2,-0.5,-0.55, -0.5,1)
			#cat(paste("***ref.index int:", initsLM$ref.index, ";length of :", length(initsLM$inits),"\n",sep=""));	
			fitShift<-fitShifts5plSingle(ydata,x, c(a,d,gainAdjust.fit$par[1:3]), data.aggregated, 
					debug=debug,
					ref.index=initsLM$ref.index, initsLM=initsLM$inits
			);
		#}
			
			if(debug){
			gainAdjust.plot(fitShift,ydata, x,ylog=data.ylog , ref.index=initsLM$ref.index,
				model="fix.all",aggregated=data.aggregated,
				a=a, d=d, xmid=gainAdjust.fit$par[1], scal=gainAdjust.fit$par[2], g=gainAdjust.fit$par[3]
				#,filename=paste("fit2_blk_",j,".jpg", sep="")
				)#a=1.2440093153044
			
				
			x.aligned<-gainAdjust.alignData(fitShift,x, model="fix.all", ref.index=initsLM$ref.index)##<---for 4Pl, fixing highest end only

			input.aligned<-gainAdjust.prepareInput(ydata,x.aligned)
			yinput.aligned<-input.aligned$y
			xinput.aligned<-input.aligned$x		
			jpeg(filename=paste(arrayName,"fit2_dot_blk_",j,"%d.jpg", sep=""))		
			plot((xinput.aligned), (yinput.aligned), type="p", main="5-p logistic fitting of intensity vs. PMT gain",
				xlab="PMT gain (log)",ylab="intensity (log)",log="xy"
				)
			#plot(log(xinput.aligned), (yinput.aligned), type="p")
			xsim<-seq(4,7,by=0.01)
			
			lines(exp(xsim), (f5pl(c(a,d, gainAdjust.fit$par[1:3]),xsim)),lty=1,col=3,lwd=2)
			dev.off();
			}			
		#done for second try
			
		#now do the last fitting for parameters, mainly for parameter a
		#-------->now fit 4 parameters on aligned data again to allow a to vary
			ydata_ag<-gainAdjust.aggregate(ydata, mode=aggregated.mode)
			aggregated_3rd<-T
			x.aligned<-gainAdjust.alignData(fitShift,x, aggregated=T, model="fix.all", ref.index=initsLM$ref.index)##<---for 4Pl, fixing highest end only

			input.aligned<-gainAdjust.prepareInput(ydata_ag,x.aligned)
			yinput.aligned<-input.aligned$y
			xinput.aligned<-input.aligned$x		
			parStartFpl<-c(a,
							gainAdjust.fit$par[1],gainAdjust.fit$par[2],gainAdjust.fit$par[3])
			cat("\tfinal fitting for parameters.....\n")
			gainAdjust.fitFinalP<-nls.lm(par=parStartFpl, 
					fn=gainAdjust.fnRes4pFpl,
					#fn=gainAdjust.fnRes3pFpl,
					#y=log(yinput.aligned), 
					y=yinput.aligned, 
					x=log(xinput.aligned),#xlen=xlen, aggregated=aggregated_3rd,
					ylog=ylog, model.weight="power"#"uniform"#"log"#"sqrt"
					,order=1,d=d,#a=0.1,
					control = nls.lm.control(nprint=1, maxit=100)
			)
			if(gainAdjust.fitFinalP$par[1]<0)
			{
				#gainAdjust.fitFinalP$par[1:4]<-c(a, gainAdjust.fit$par[1:3])
				gainAdjust.fitFinalP$par[1]<-min(ydata)
			}
			if(debug){
				gainAdjust.plot(fitShift,ydata_ag, x,ylog=data.ylog , ref.index=initsLM$ref.index,
						model="fix.all",
						a=gainAdjust.fitFinalP$par[1], 
						,d=d, 
						xmid=gainAdjust.fitFinalP$par[2], scal=gainAdjust.fitFinalP$par[3], g=gainAdjust.fitFinalP$par[4],
						aggregated=T
						#,filename=paste("fit3_blk_",j,".jpg")
				)#a=1.2440093153044
				
				jpeg(filename=paste(arrayName,"fit3_dot_blk_",j,"%d.jpeg",sep=""))
				plot((xinput.aligned), (yinput.aligned), type="p", main="5-p logistic fitting of intensity vs. PMT gain",
					xlab="PMT gain (log)",ylab="intensity (log)",log="xy"
					)
				#plot(log(xinput.aligned), (yinput.aligned), type="p")
				xsim<-seq(4,7,by=0.01)
				
				lines(exp(xsim), (f5pl(c(gainAdjust.fitFinalP$par[1],d, gainAdjust.fitFinalP$par[2:4]),xsim)),lty=1,col=3,lwd=2)
				dev.off();
			}
			cat("\tDONE!\n");
			#shift<-
			#fitRes<-list("5PL"=c(gainAdjust.fitFinalP$par[1],d, gainAdjust.fitFinalP$par[2:4]), 
			#	"shift"=gainAdjust.insertAt(fitShift$par,0, initsLM$ref.index))
			#fr[[j]]<-fitRes;
			
			#from this point one, we start doing the "normalization" to pick the good ones for each array
			#first do control
			if(fit.mode=="all")
			{
				ydata<-y$C[idx.prC,]
				initsLM<-gainAdjust.prepareInitsLM(ydata,x, aggregated=data.aggregated)
			}
			#get shift
			#initsLM<-gainAdjust.prepareInitsLM(ydata,x, aggregated=data.aggregated)
			if(j!=1)  ###this is not necessary, since inside prepareInitsLM it has been set to be -1 when using outside reference line
			{			###but keep it there anyway to make sure.
				initsLM$ref.index<- -1
			}
			fitShift<-fitShifts5plSingle(ydata,x, c(gainAdjust.fitFinalP$par[1],d,gainAdjust.fitFinalP$par[2:4]), 
											data.aggregated, ref.index=initsLM$ref.index, initsLM=initsLM$inits);
			if(j==1)  #fitShift do not contains the ref.index shift (which is zero anyway), so when we use
			{			#the inside referene (at the first round of estimate) we have to insert it in the shift
						#because we will use it to generate all the data in "adjust.Matrix" function
				fitShift$par<-gainAdjust.insertAt(fitShift$par, 0, initsLM$ref.index)
			}
			cshifts[[j]]<-fitShift$par #gainAdjust.insertAt(fitShift$par, 0, initsLM$ref.index)
			#generate output
			temp_gain<-NULL
			if(!missing(PMT.gain)&&!is.null(PMT.gain)){
				temp_gain<-PMT.gain$C[idx.prC]
			}
			#cat("1\n")
			ydataC.adj<-adjust.Matrix(ydata=ydata,x=x,par.5pl=c(gainAdjust.fitFinalP$par[1],d,gainAdjust.fitFinalP$par[2:4]),fitShift=fitShift,PMT.gain=temp_gain, F);
			#cat("2\n")
			fr$C[idx.prC]<-ydataC.adj$E;
			fr$C.adj[idx.prC]<-ydataC.adj$E.adj;
			fr$C.gain[idx.prC]<-ydataC.adj$gain;
			#cat("3\n")
			
			#----------------
			#Now do fitting for shifts for target expression only
			idx.prE<-which(is.element(y$gene[,"Block"],c(bg:ed)))
			ydata<-(y$E[idx.prE,])
			#cat("4\n")
			initsLM<-gainAdjust.prepareInitsLM(ydata,x, aggregated=data.aggregated)
			if(fit.mode=="control")
			{
				#here we are doing target/feature ones, but we are doing control fitting, so we need to 
				#make the ref.index in this case to be external, not internal
				
				#initsLM$inits<-gainAdjust.insertAt(initsLM$inits,0,initsLM$ref.index)
				initsLM$ref.index<- -1
			}
			temp_gain<-NULL
			fitShift<-fitShifts5plSingle(ydata,x, c(gainAdjust.fitFinalP$par[1],d,gainAdjust.fitFinalP$par[2:4]), data.aggregated, ref.index=initsLM$ref.index, initsLM=initsLM$inits);
			
			#fitShift$par<-gainAdjust.insertAt(fitShift$par, 0, initsLM$ref.index)
			
			shifts[[j]]<-fitShift$par #gainAdjust.insertAt(fitShift$par, 0, initsLM$ref.index)
			if(!missing(PMT.gain)&&!is.null(PMT.gain)){
				temp_gain<-PMT.gain$E[idx.prE]
			}
			ydataE<-adjust.Matrix(ydata=ydata,x=x,par.5pl=c(gainAdjust.fitFinalP$par[1],d,gainAdjust.fitFinalP$par[2:4]),fitShift=fitShift,PMT.gain=temp_gain,F);
			fr$E[idx.prE]<-ydataE$E;
			fr$E.adj[idx.prE]<-ydataE$E.adj;
			fr$E.gain[idx.prE]<-ydataE$gain;
			fit[[j]]<-c(a=gainAdjust.fitFinalP$par[1],d=d, xmid=gainAdjust.fitFinalP$par[2],scal=gainAdjust.fitFinalP$par[3],g=gainAdjust.fitFinalP$par[4])
			
			#shifts[[j]]<-gainAdjust.insertAt(fitShift$par,0,  )
			#cat("show fit paraemter.....:>>>>\n")
			#print(fit[[j]])
			#cat("\n")
		} #end of for loop for analysis by block
	
	#prepare the output
		fr$fit<-fit
		fr$shifts<-shifts
		fr$cshifts<-cshifts
		#cat("--------show all fit paraemter for one array.....:>>>>\n")
		#	print(fr$fit)
		#	cat("\n")
		fr
	}#end of the function
	
	#'@title PMT gain adjust data
	#'@description use fitted 5-parameter logistic function to adjust the 
	#'	elist protoArray data
	#'@details In this function, based on the fitted 5-parameter logistic function (5pl),
	#'	the elist protoarray data is adjust for each feature/control so as to select the 
	#'	signals with the PMT gain closest to the middle point of curve 
	#'	(xmid) and are aggregated across the duplicates by predict  
	#'	the selected PMT gain to the fitted 5pl. 
	#'@seealso \code{\link{gainAdjust.prepareInitsLM}}
	#'
	#'@param ydata data matrix for either elist expression field, "E" or control signal
	#'		field, "C". This is the raw data to be adjusted.
	#'@param x numeric vector for the PMT gain settings for generating the data
	#'@param par.5pl the non-linear fitting of 5pl. It contains the parameter fitted
	#'		using the data lines.
	#'@param fitShift numeric vector containing the fitting shifts for each feature line.
	#'		it contains the shift for the reference lines if the fit is done by an internal
	#'		reference, of which the shift is zero.
	#'
	#	fitShift including all the shifts for all data points for each individual Array data E or C
	#	This function includes the one with index of initsLM$ref.index.
	#'@param PMT.gain a numeric vector it has the identical order as to that of feature
	#'	in ydata. It could be missing or null, which means it has not been determined yet.
	#'	In this case, we need to determine the best PMT gain for each feature and pass it
	#'	as output.
	#'
	#'@param data.aggregated logic used to indicate whether the data haven't aggregated
	#'	
	#'@return a data list contains five fields: E field, the unmodified signal data at the
	#'	best PMT gain ( closest to the xmid ); E.adj field, the predicated signal value by 5pl at
	#'	the best PMT gain for each feature. This could also be used to aggregated 
	#'	the duplicated data and will be the same for the duplicated feature data. ; PMT.gain,
	#'	the numeric vector recording the best PMT gain for data value. If PMT.gain is set
	#'	as input, then we simply use these values and pass it as output too. If otherwise,
	#'	we need to generate PMT.gain based on the fitted 5pl function and pass it as output.
	#'
	#'@export
	adjust.Matrix<-function(ydata, x, par.5pl,  fitShift, PMT.gain, data.aggregated=F )
	{
		if(!data.aggregated)
		{
			fitShift$par<-rep(fitShift$par, rep(2, length(fitShift$par)))
		}
		##check for the data integrity
		if(missing(PMT.gain)||is.null(PMT.gain)){#this is the case we need to compare to get optimal pmt for each gene
			PMT.gain<-fitShift$par
			#now go through the list to pick the  best one
			xmid<-par.5pl[3]
			for(i in 1:length(fitShift$par))
			{
				PMT.gain[i]<-x[which.min(abs(log(x)+fitShift$par[i]-xmid))]
			}
		}
		#now for sure, we have PMT.gain reference setting, get 
		y.out<-list("E"=c(),"E.adj"=c(), "gain"=PMT.gain)
		y.out$E.adj<-f5pl(par.5pl,(log(PMT.gain)+fitShift$par))
		y.out$E<-fitShift$par
		for(j in 1:dim(ydata)[1])
		{
			y.out$E[j]<-ydata[j,which(x==PMT.gain[j])]
		}
		y.out
	}
	
	#function to insert value at "index" into an array
	#'@export
	gainAdjust.insertAt<-function(vec, value=0, index=1)
	{
		if(missing(vec))
		{
			stop("missing the input vector, please specify")
		}
		if(class(vec)!="numeric")
		{
			stop("the input vector is not in correct format. please specify")
		}
		
		arr<-c(0, vec)
		if(index==1)
		{
			arr[1]<-value
		}else if(index==length(vec)+1)
		{
			arr[1:length(vec)]<-vec
			arry[length(vec)+1]<-value
		}else  #in this case, index is not on th boundary
		{
			arr[1:(index-1)]<-vec[1:(index-1)]
			arr[index]<-value
			arr[(index+1):length(arr)]<-vec[index:length(vec)]
		}
		arr
	}
	
	###obsoleted Not used  <----- 
	###
	##assuming has done background correction
	##input, 
	###	x, is the PMT gain setting vector/matrix assuming all the arrays having the identical pmt vector
	##Obsoleted NOt used 
	#  # '@export
	gainAdjust.fitArrayByBlock<-function(object, x, ylong=F, aggregated=F, debug=F )
	{
		if(missing(object))
		{
			stop("missing object, please specify....")
		}
		if(class(object)!="EListRaw")
		{
			stop("object is not of correct format. please specify......")
		}
	
		#now getting data ready to do the fitting......
		blkNum<-max(object$gene[,"Block"])
		arrys<-unique(object$target[,"Array"])
		arryNm<-length(arrys)
		if(class(x)=="numeric")
		{
			x<-matrix(rep(x, arryNm),nrow=arryNm, ncol=length(x),byrow=T)
			#x<-
		}
		#check for the dim
		if(dim(x)[1]<arryNm)
		{
			stop("x is not in correct form, please specify!")
		}
		cat("Start doing the fitting ......\n")
		for(i in 1:blnNum) #by blocks first
		{
			cat("block ", i,"/",blnNum, "...\n")
			for(j in 1:arrys) #by array 
			{
				#now take out the array data from the object
				aNm<-arrys[j]
				idx.aNm<-which(object$target[,"Array"]==aNm)
				#now we get the index of the one array j, now need to do fitting by block
				
				#for the array, we need to get the protein data by block
				idx.prE<-which(is.element(object$gene[,"Block"],c((i+4):(i+7)))) #which(object$gene[,"Block"]==i|object$gene[,"Block"]==i+1|object$gene[,"Block"]==i+2|object$gene[,"Block"]==i+3)
				idx.prC<-which(is.element(object$cgene[,"Block"],c((i+4):(i+7))))  #which(object$cgene[,"Block"]==i|object$cgene[,"Block"]==i+1|object$cgene[,"Block"]==i+2|object$cgene[,"Block"]==i+3)
				
				#now get idx, need to get the data into matrix
				y<-rbind(object$E[idx.prE, idx.aNm], object$C[idx.prC,idx.aNm])
				x.arry<-x[j,]
				
				#now doing fitting
				f<-gainAdjust.fit5PL(y, x.arry,F, aggregated=F, debug=T)
				
				plot(c(250,650),c(0.05, 65553),log="xy",type="n")
				for(k in 1:length(y[,1]))
				{
					lines(x.arry,(y[k,]), col=k)
				}
			}
		}
	
	}
	
#'@title PMT gain adjust of protoarray data with a 5 parameter logistic function 	
#'@description 
#' this is the outermost driving function to do the fitting and adjustment on EListRaw data of protoArray
#'@details please see details of \code{\link{gainAdjust.fit5plArray}}
#'@return data list contains the following fields
#' \itemize{
#'		\item{E, the raw feature data at the best PMT gain determined by 5pl}
#'		\item{C, the raw control data at the best PMT gain determined by 5pl}
#'		\item{E.adj, the adjusted/predicted feature data at the best PMT gain determined by 5pl}
#'		\item{C.adj, the adjusted/predicted control data at the best PMT gain determined by 5pl}
#'		\item{C.PMT, the the best PMT gain for control data determined by 5pl}
#'		\item{E.PMT, the the best PMT gain for feature data determined by 5pl}			
#'		\item{fit, the obtained 5pl fitting }
#'		\item{shift, the shifts for feature data determined by 5pl}
#'		\item{cshift, the shifts for control data determined by 5pl}
#'}
#'@seealso \code{\link{gainAdjust.fit5plArray}}
#'@export
PMT.adjust<-function(data, #raw date after reading the GPR files
			x, #vector holding the PMT gains for reading the chips, the arrangement of data in PMT reading is identical to the element order in x.
			data.ylog=F, data.aggregated=F,#description of data, usually, don't need to change, for advance user only
			fit.mode="control", ##fit the 5pl with control data only, to speed up 
			fit.aggregated=F, #fit the 5pl with aggregated data, for advanced data
			aggregated.mode="geometric", #controlling the way to do aggregation, geometric is preferred
			block.size=12, #size of blocks to do fitting, the more the better??. but too slow, 12 is appropriate
			threshold=60,
			a=1.1, d=65535, #change only when you know what you do.
			debug=T  #showing the debugging information
			)
{
	#check data
	if(missing(data)||missing(x))
	{
		stop("please specify the data/x.")
	}
	if(class(data)!="EListRaw")
	{
		stop("the data object is not in correct format!!")
	}
	if(class(x)!="numeric")
	{
		stop("the x is not in correct format!!")
	}
	if(floor(dim(data$E)[2]/length(x))*length(x)!=dim(data$E)[2])
	{
		stop("data and x are not correctly set up, please check")
	}
	#now start picking the data 
	blkNum<-max(data$gene[,"Block"])
	arrys<-unique(data$target[,"Array"])
	arryNm<-length(arrys)

	#i<-j<-1
	#cat("block ", i,"/",blkNum, "...\n")
	y<-list("C"=NULL, "E"=NULL, "cgene"=NULL, "gene"=NULL)
	yAdj<-list("E"=matrix(0,ncol=arryNm, nrow=dim(data$E)[1],byrow=T), 
				"C"=matrix(0,ncol=arryNm, nrow=dim(data$C)[1],byrow=T), 
				"E.adj"=matrix(0,ncol=arryNm, nrow=dim(data$E)[1],byrow=T), 
				"C.adj"=matrix(0,ncol=arryNm, nrow=dim(data$C)[1],byrow=T),
				"E.PMT"=NULL, "C.PMT"=NULL,
				"fit"=vector(mode="list", length=arryNm),
				"shifts"=vector(mode="list", length=arryNm),#shifts for 5pl fitting, for targets/features,
				"cshifts"=vector(mode="list", length=arryNm)#shifts for 5pl fitting, for control proteins
			)

	PMT.gain<-NULL;
	
	for( j in 1:arryNm)
	{
		#now take out the array data from the object
		aNm<-arrys[j]
	
		idx.aNm<-which(data$target[,"Array"]==aNm)
		#now we get the index of the one array j, now need to do fitting by block
		
		#now get idx, need to get the data into matrix
		
		y$C<-data$C[,idx.aNm] ; #rbind(object$E[idx.prE, idx.aNm], object$C[idx.prC,idx.aNm])
		y$E<-data$E[,idx.aNm];
		y$cgene<-data$cgene;
		y$gene<-data$gene;
		
		#if(j==1)
		#{
		#	
		#}
		rs<-gainAdjust.fit5plArray(y, x, data.ylog=data.ylog, data.aggregated=data.aggregated,
		fit.mode=fit.mode,
		fit.aggregated=fit.aggregated, aggregated.mode=aggregated.mode, order=1,
		block.size=block.size,a=a, d=d,threshold=threshold,
		PMT.gain=PMT.gain,arrayName=aNm,
		debug=debug
		)
		#cat(paste("\n===done for one array",arrys[j],"\n",sep=""));
		#need to get PMT.gain using the first one as reference
		PMT.gain<-list("E"=rs$E.gain, "C"=rs$C.gain);
		
		#at end of each run, need to add the result to the output list
		yAdj$E[,j]<-rs$E	
		yAdj$C[,j]<-rs$C
		#colnames(yAdj$C)[j]<-as.character(arrys[j])
		yAdj$E.adj[,j]<-rs$E.adj
		#colnames(yAdj$E.adj)[j]<-as.character(arrys[j])
		yAdj$C.adj[,j]<-rs$C.adj
		#colnames(yAdj$C.adj)[j]<-as.character(arrys[j])
		yAdj$fit[[j]]<-rs$fit
		yAdj$shifts[[j]]<-rs$shifts
		yAdj$cshifts[[j]]<-rs$cshifts
	}
	colnames(yAdj$E)<-as.character(arrys)
	colnames(yAdj$C)<-as.character(arrys)
	colnames(yAdj$E.adj)<-as.character(arrys)
	colnames(yAdj$C.adj)<-as.character(arrys)
	yAdj$E.PMT<-PMT.gain$E
	yAdj$C.PMT<-PMT.gain$C
	yAdj;
}

####code to plot/visualize the fitted lines by block
#
#data, matrix, the expression/signal data matrix for one single array, for example E or C table
#		if missing then we only plot the fitted line without data
#		The data is arranged in such a way that the row is for different genes/proteins in an order 
#		identical to the genes table and the column is for different PMT gains for the same array.
#pmts, vector, listing the different PMT gains for the data
#genes, data.frame the gene table used to disecting data into blocks to fit for the 5pl lines
#
#block.size, scalar number, size for the data to fit the line.
#
#fit, the list for the 5pl parameters. each fit is a array of 5 parameters
#		the total number of fit parameters is determined by the block.size
#
#'@export
gainAdjust.plotCompareFit<-function(data,pmts, genes,block.size=12, fit, shifts, data.aggregated=F)
{
	if(missing(fit)||length(fit)==0)
	{
		stop("please specify the fitting parameters")
	}
	#first determine some parameters first
	xsim<-seq(4,7,by=0.1)
	ylim<-c(fit[[1]][1],fit[[1]][2])
	
	#
	#manipulate the data and make it 
	#
	
	if(!missing(data))
	{##we are plotting the data lines
		if(missing(pmts))
		{
			stop("please specif the PMTS gains for the data")
		}
		if(missing(genes))
		{
			stop("Please specify the genes table for the data")
		}
		#if(missing(shifts))
		#{	
		#	stop("please specify the shifts for the data lines")
		#}

		xsim<-seq(min(min(log(pmts)),4),max(max(log(pmts)),7),by=0.1)
		ylim<-c(min(fit[[1]][1],min(data)),max(fit[[1]][2], max(data)))
		
	}
	plot(c(min(xsim),max(xsim)), ylim, type="n", main="5-p logistic fitting of intensity vs. PMT gain",
					xlab="PMT gain (log)",ylab="intensity (log)",log="y"
					)
	if(!missing(data))
	{
		#parse the data and make them read
		blkNum<-max(genes[,"Block"])
		for(j in 1:ceiling(blkNum/block.size)){
			cat("plotting by blocks: ",j, "/",blkNum/block.size, "...\n")
			#for the array, we need to get the protein data by block, 
			bg<-(j-1)*block.size+1
			ed<-(j-1)*block.size+block.size
			#cat(paste("bg:",bg, "; ed:", ed,"\n", sep=""));
			if(ed>blkNum)
			{
				ed<-blkNum
			}
			
			idx.prC<-which(is.element(genes[,"Block"],c(bg:ed)))
			ydata<-data[idx.prC,]
			shifts_blk<-shifts[[j]]
			if(!missing(shifts))
			{
				#if(!data.aggregated)
				#{
				#	tms<-rep(2,length(shifts))
				#	shifts<-rep(shifts,tms)
				#}
				for(i in 1:length(shifts_blk))
				{
					if(data.aggregated){
						lines(log(pmts)+shifts_blk[i],ydata[i,],col=j, lty=j)
					}else{
						lines(log(pmts)+shifts_blk[i],ydata[i*2,],col=j, lty=j)
						lines(log(pmts)+shifts_blk[i],ydata[i*2-1,],col=j, lty=j)
					}
					
				}
			}
		}
	}
	
	#plot(log(xinput.aligned), (yinput.aligned), type="p")
	#xsim<-seq(4,7,by=0.01)
	#now start doing the	
	for(j in 1:ceiling(blkNum/block.size)){
		lines((xsim), (f5pl(fit[[j]],xsim)),lty=j,col=j,lwd=2)
	}
	
	#done
}				

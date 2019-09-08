##=======gainAdjust_PowerModel.R==========
#
## This is the module to do the gainAdjust using the newer PMT exponent model. 
##The job is based on the PMT gain model (see the explanations in the notes) Tom's note and my notes. 
##We start from readings under multiple gain settings, "summarize" the reading in order
## to estimate the true light influx
## 
## The model is based on a baseline added exponential model (log-transformed exponential with baseline).
## 
## We assume all the readings should be well described a 
## common parameterized function. The different
## protein concentration as well as the different affinity
## together determine/shift the initial light influx (scaled) and then
## the other parameters, baseline, saturated intensity and exponent
## identical to all the spots on the same array. 
## . we will fit k light influx (k is the number of spots, count repeats as one spot)
## plus 2 other parameters baseline and exponent.
## Another thing to note is that we assume the maximum intensity is 
## 65535, but we will fit the baseline/background.
## 
##	B=0.01, d=65535
## B is almost zero, because we will do background correction before doing the fitting.
## another thing to mention: the form of the  function we assume
## here is the one assume a log transformed x input, since it takes an exponetial
## in the formula. We will all try to using the regular scale.
## 
## For y, we will take the log, since we assume the residual is logError.
## actully, we don't have to do this, since we can use the regular scale, but we fit
## the nonlinear using the weight to control the residual error. 
## 
## 	We will try both anyway
##
###Note: logistic function 5pl gainadjust is done in a different module "gainAdjust.R"
#   Model: log(y)<-log(B+fi(v)^delta)+eps
#for the 
###model 
##  log(y)<-log(B+(v*phi)^delta)+eps
### y - log intensity
### B - background
### phi - (fi)^(1/delta)
### delta - exponent
### eps  - residual error
### 
#####

### model log(x)
###
##  log(y)<-log(B+exp((log(v)+log(phi))*delta))+eps
### see above
#####
##=======created 1/25/2019 
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
pmt_saturated<-65535;
pmt_baseline<-log(0.1);
minimum.phi<- -6.5;#-7 previously #smallest possible phi value. it is not possible to have a value estimated smaller than this. it is possibly due to errors
#B<-pmt_lowest_reading;
#G<-(pmt_saturated_reading);

#'@include residualFunc.R
	
	#'@title prepare initial values  
	#'@description generate the initial values for fitting scaled light influx with 
	#'	a PMT gain power model.
	#'@details this is a functoins to prepare the initials for the exponent of the PMT gain power model. 
	#'	The assumption is that the intensity of a feature in protein microarray
	#'	(ProtoArray, invitrogen) acquired with different PMT gains follow 
	#'	a power function  as \cr
	#'	\ifelse{html}{\out{<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-MML-AM_CHTML"></script>
	#'						$$y = B + {(phi*v)}^{delta} $$
	#'				 	}}{\deqn{y=B+ {(phi*v)}^{delta}}{ASCII}}
	#'	\cr 
	#'
	#'	where a is the upper bound; d is the lower bound);
	#'	or alternatively \cr
	#'  \ifelse{html}{\out{<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-MML-AM_CHTML"></script>
	#'						$$y = B +  phi*{v}^{delta} $$
	#'				 	}}{\deqn{y=B+  phi*{v}^{delta}}{ASCII}}
	#'	\cr 
	#'	or
	#'   \ifelse{html}{\out{<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-MML-AM_CHTML"></script>
	#'						$$log(y) = log(B + exp(xmid*scal+ log(v)*scal)) $$
	#'				 	}}{\deqn{log(y)=log(B+ exp(xmid*delta+ log(v)*scal))}{ASCII}}
	#'	\cr 
	#'	in this package, we use the first one in our implement (avoiding the overflow)\cr
	#'	On the same array, or across different array, we have identical delta and B, but 
	#'	different phi for each spot/feature.
	#' \cr
	#'	The algorithm takes a local log-linear model to generate the initial values like this,
	#'	
	#'	\itemize{
	#'	\item {1. find the largest two intensity in each serier, the second largest one of which is smaller than 
	#'		saturated intensity\cr}
	#'		{
	#'		if all intensities are no smaller than the saturated intensity, we need to label this as
	#'			two larget for estimatoin (??? what to do???).
	#'		}
	#'	\item {2. determine the initial k values using log linear assuming there is no baseline\cr}
	#'		{
	#'		slope is the delta, and fi is exp(intercept/slope) 
	#'		
	#'		}
	#'	\item {3.the ouput is a list with a vector of all the initials, and a list of index pointing
	#'		to ones that are all saturated(??)\cr}
	#'		{
	#'		the saturated ones are no information to determine the accurate fi, we need to set it to some value
	#'		Hopefully, there will not be many cases like this.
	#'		}
	#'}
		
	
	#'@param y dataframe this is raw data and has NOT been prepared. 
	#'		it has row as genes and col as pmts
	#'@param x numeric vector It contains the PMT settings for each array. 
	#'@param data.aggregated Parameter indicates whether the data has been 
	#'		aggregated, mean the two repeats has been averaged. 
	#'		By default it is not. In this case, we will generate 
	#'		the identical shift, k, for the two repeats.
	#'@param data.log indicating whether the data has been log-transformed. 
	#'@return a data list contain two things. First is vector of all initial values
	#'	it has the identical length to the spots (will need to aggregate repeats.
	#'	The second is index pointing to the values where we all has saturated intensity for 
	#'	all PMT gain setting.
	#'		Note: 1) the prepared inits vector always is in a format of aggregated data
	#'		; 2)the initial value of phi is not the final fomat of phi, which is (phi)^delta.
	#'		3)also when we want to use these phi values as the initials, we need to 
	#'		log-transformed to feed the fitting function.
	#'@seealso [gainAdjust_fit_Pbp_phi()][gainAdjust_fit_Pbp_array()][gainAdjust_fit_Pbp_data]
	#'		[f_PowerBaselinePMT()]
	#'@export
	gainAdjust_prepareInits_Pbp<-function(y, x, aggregated=TRUE, data.log=F)
	{
		if(missing(y))
		{
			stop("ERROR:missing input y, please specify.........");
		}
		if(missing(x))
		{
			stop("ERROR:missing input x, please specify.........");
		}
		if(class(y)!="matrix")
		{
			stop("the input y is expected to be matrix, please check.......");
		}
		if(class(x)!="numeric")
		{
			stop("the input y is expected to be numeric vector, please check.......");
		}
		#get dimensions
		rlen<-dim(y)[1]
		clen<-dim(y)[2]
		yTemp<-y;
		xTemp<-x;
		if(!data.log)
		{
		 yTemp<-log(y)#y[seq(1,rlen,by=2),];
		 xTemp<-log(x)
		}
		
		
		####need to check for the aggregated status, in order to "arrange" data
		if(!aggregated)#need to summarize the data
		{	
			yTemp.agg<-yTemp[seq(1,rlen,by=2),]
			if(rlen%%2!=0)
			{
				stop("the number of array rows are not even, can not do aggregation! Stop!! ");
			}
			#for(i in seq(1, rlen, by=2))
			#{
			#	yTemp[(i+1)/2,]=(y[i,]+y[i+1,])/2
			#}
			in.y<-seq(1,rlen,by=2)
			yTemp.agg[(in.y+1)/2,]=(yTemp[in.y,]+yTemp[in.y+1,])/2
			#if(!missing(ref.line))
			#{
			#	ref.line$y[1,]=(ref.line$y[1,]+ref.line$y[2,])/2
			#}
			yTemp<-yTemp.agg;
		}
		
		rlen<-dim(yTemp)[1]
		#ref<-yTemp[1,clen]
		
		#sort x in case X is not in order
		idx.srtX<-order(xTemp)  #increasing order
		
		inits<-rep(0,rlen)
		Y<-(yTemp[,idx.srtX])
		if(class(Y)!="matrix")
		{
			if(class(Y)!="numeric")
			{
				stop("the input y is expected to be matrix, but something wrong, please check........");
			}
			Y<-matrix(Y, nrow=1, byrow=T);
		}
		X<-xTemp[idx.srtX]
		
		idx.saturated<-c(); #will be used to hold the indexes of lines that are saturated points.
				#which.min(abs(pmt_saturated_reading/2-Y[,clen]))##here yTemp is the sum of the duplicated data, so we choose the reference to be the one closest to MaxF/2
		#idxs.arr<-yTemp[-1])
		#ref<-yTemp[idx.ref,idx.srtX]
		#if(!missing(ref.line))
		#{
		#	idx.ref<- -1
		#	ref<-ref.line$y[1,idx.srtX]
		#}
#		cutoff<-(pmt_saturated*0.1); #used to define a good unsatured point
		delta<-rep(0,rlen);
		phi<-rep(0,rlen);
		for(i in c(1:rlen))
		{
			#do the job to calculate the initial values based on log-linear ks
				#first find the largest Y pair, which are smaller than saturated values. 
			#want to find a Yi that is smaller than ref_max
			Yi<-Y[i,]
			
#			j<-1;
#			Y1<-Yi[length(Yi)];
#			X1<-0;X2<-0;
#			Y2<-Yi[length(Yi)-1];
#			for(j in (length(Yi)):2)
#			{
#				Y1<-Yi[j];
#				Y2<-Yi[j-1];
#				X1<-X[j];
#				X2<-X[j-1];
#				if((Y1<log(pmt_saturated-cutoff)&&Y2<log(pmt_saturated-cutoff))&&Y1>Y2)
#				{
#					break;
#				}
#			}
#			if(Y2>log(pmt_saturated-cutoff)) #in this case, we are having all saturated points
#			{
#				cat("find one.........", "i:",i,"Y1:",Y1,";Y2:",Y2,"\n");
#				idx.saturated<-c(idx.saturated,i);
#			}

			#now start doing matrix regression
			b<-regressionByMatrix(Yi, X)
			#whether the line is saturated, we still do it.
#			slope<-(Y1-Y2)/(X1-X2)
			slope<-b[2];
			
#			intercept<-((Y1/X1)-(Y2/X2))/(1/X1-1/X2)
					#Xi_r_pred<-(Yi_j-b)/a
			intercept<-b[1]
			delta[i]<-slope;
			phi[i]<-exp(intercept/slope);
				
		}##for loop for each line/series
		phi[phi>0.009]<-median(phi); ###this is one is get rid of the one that is not right due to 2 point estimation
		phi[phi<=0.0025]<-median(phi);
		phi[phi<=0.0025]<-0.0025;
		phi[phi>0.009]<-0.008;
		#list("inits.delta"=median(delta), "inits.phi"=phi,"index.saturated"=idx.saturated);
		list("inits.delta"=(sort(delta)[floor(length(delta)*0.9)]), "inits.phi"=phi,"index.saturated"=idx.saturated);
	}
	
	#'@title plot the fitted and raw data of the PMT power model with baseline.
	#'@description plot gene signals versus PMT voltage and the fitted data. 
	#'@details the fitted data are the results from a fitted line with the PMTs
	#'	power model with baseline. No shifting the lines.
	# #  ##'@param B a scalar denoting the lower background
	#'@param G a scalar denoting the saturated intensity, 2^16-1 for 16bit encoding.
	#'@param phi the estimating scaled light influ
	#'@param delta the model exponent 
	#'@param log string indicating whether to do log scale while plotting
	#'@param ref.index the index of the refer line. All other lines will be "shifted" toward it.
	#'@param filename for the output plot. if not provided, will be output to the stdout 
	#'@seealso [f_PowerBaselinePMT()]
	#'@export  
	#<--------------to do 
	gainAdjust_plot_Pbp_raw<-function(LMfit,y, x, data.aggregated=FALSE,
			#model="fix.both",
			#B,
			G=pmt_saturated, 
			#delta, phi, ref.index=1,
			log=c("xy","","x","y"),
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

		if(missing(y))
			{
				stop("the input y is missing, please specify!!!")
			}
		if(missing(x))
			{
				stop("the input x is missing, please specify!!!")
			}
		if(!data.aggregated)
		{
			y<-gainAdjust.aggregate(y, mode="geometric");
		}
		
#		switch(model,
#				"fix.both"={	
#						#this has been taken care of as the default case in above
#						xmid<-LMfit$par[1]
#						scal<-LMfit$par[2]
#						g<-LMfit$par[3]
#						#do it again here, because we can not leave it empty
#						k_index<-4
#					},
#				"fix.high"={
#						#be careful here, 4PL is not the real 4PL, it is 5PL but fixed d, maximum level
#						a<-LMfit$par[1]
#						xmid<-LMfit$par[2]
#						scal<-LMfit$par[3]
#						g<-LMfit$par[4]
#						k_index<-5
#					},
#				"fix.none"={
#						#all five parameter
#						a<-LMfit$par[1]
#						d<-LMfit$par[2]
#						xmid<-LMfit$par[3]
#						scal<-LMfit$par[4]
#						g<-LMfit$par[5]
#						k_index<-6
#				},
#				"fix.all"={
#					if(missing(xmid)|missing(scal)|missing(g))
#					{
#						stop("please specify parameters (xmid/scal/g)!!!")
#					}
#					k_index<-1
#					
#				},
#				stop("unkown model specified")
#			)
		#now doing the plot
		if(!is.null(filename)&&(filename!=""))
		{
			jpeg(filename=filename,quality=98, res=300, width=1600, height=1500)
		}
		#op<-par(mfrow=c(2,1), mar = c(3,5,2,1))
		#need to figure out the range of y and x
		#for y, we know it is (0.1, 65535)
		
		#for x, we need to add ks
		xmin<-min(x)
		xmax<-max(x)
		
		ymin<-min(y)
		ymax<-max(y)
		
		#take out the parameters
		B<-LMfit$par[1];
		delta<-LMfit$par[2];
		phi<-LMfit$par[-c(1,2)];
		if(log!="")
		{
			log<-match.arg(log);
		}
		plot(c(xmin*0.95,xmax*1.05), c(ymin*0.95,ymax*1.05), bg = "black", #cex = 0.5, 
				main="PMT Power Model Fitting", type="n",log=log
				, xlab="PMT voltage(v)", ylab="Intensity"
				)
		#plot(c(xmin,xmax+0.3), c((0.1),(66535)+5000), bg = "black", cex = 0.5, main="data", type="n")
		#points(log(x), log(y[1,]), col="grey", lty=2) 
		
		#now plot the fitted value
		xx<-seq(xmin,xmax, by=(xmax-xmin)/1000)
		
		for(i in c(1:length(phi)))
		{
			points(x, y[i,], col=i, cex=0.7)
			yy<-f_PowerBaselinePMT(c(B, G, phi[i], delta),xx);
			lines(xx,yy,col=i)
		}
		if(!is.null(filename)&&(filename!=""))
		{
			dev.off();
		}
	}
	
#'@title plot the fitted and raw data for the PMT power model with baseline 
#'@description plot the raw data and fitted line for debugging. The fitted model
#'	is the results of PMT power model with baseline non-linear model.
#'@details The parameters for the fitted line are in the results of the fitting. We shift 
#'	the raw data based on the phi parameter. All lines are shifted towards one reference line.
#'  Then the fitted line (reference fitted line) is plotted. Remember that phi is assumed to 
#'	be log-transformed. 
#'@param B a scalar denoting the lower background
#'@param G a scalar denoting the saturated intensity, 2^16-1 for 16bit encoding.
	#'@param phi the estimating scaled light influ, log-transformed.
	#'@param delta the model exponent
	#'@param type character indicating the plotting type. {\code{\link{plot}}}
	#'@param log string indicating whether to do log scale while plotting
	#'@param ref.index the index of the refer line. All other lines will be "shifted" toward it.
	#'@param filename for the output plot. if not provided, will be output to the stdout 
#'@seealso [f_PowerBaselinePMT()] [plot()]
#'@export	
gainAdjust_plot_Pbp<-function(LMfit,y, x, data.aggregated=FALSE,
			#model="fix.both",
			#B,
			G=gainscan:::pmt_saturated, 
			#delta, phi, 
			ref.index=1,  #used to indicate which line is used as the reference, all other ones is shifted towards it.
			log=c("xy","x","y",""), type=c("p","l","n","o","b","h","s","S","c"),
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

		if(missing(y))
			{
				stop("the input y is missing, please specify!!!")
			}
		if(missing(x))
			{
				stop("the input x is missing, please specify!!!")
			}
		if(!data.aggregated)
		{
			y<-gainAdjust.aggregate(y, mode="geometric");
		}
		
		#now doing the plot
		if(!is.null(filename)&&(filename!=""))
		{
			jpeg(filename=filename,quality=98, res=300, width=1600, height=1500)
		}
		
		if(ref.index > dim(y)[1])
		{
			ref.index<-1;
		}
		#take out the parameters
		B<-LMfit$par[1];
		delta<-LMfit$par[2];
		
		phi<-LMfit$par[-c(1,2)];
		if(log!="")
		{
			log<-match.arg(log);
		}
		if(ref.index!=1)
		{
			#first shifted the reference line to the beginning.
			temp<-y[ref.index,];
			y[ref.index,]<-y[1,]
			y[1,]<-temp;
			temp<-phi[ref.index]
			phi[ref.index]<-phi[1]
			phi[1]<-phi[ref.index]
		}
		
		#now doing the shifting only on xs
		phi.ref<-phi[ref.index];
		phi<-phi-phi[1]; #now phi[1] is the reference.
		x_shift<-(exp(phi))
		x_shift<-matrix(x_shift, ncol=dim(y)[2], nrow=dim(y)[1], byrow=F)
		x<-matrix(x, ncol=dim(y)[2], nrow=dim(y)[1], byrow=T); 
		
		x_shift<-x*x_shift
		
		#for x, we need to add ks
		xmin<-min(x_shift)
		xmax<-max(x_shift)
		
		ymin<-min(y)
		ymax<-max(y)
		
		plot(c(xmin*0.95,xmax*1.05), c(ymin*0.95,ymax*1.05), #cex = 0.1, #bg = "black",  
			main="PMT Power Model Fitting", type="n",log=log
			, xlab="PMT voltage(v)", ylab="Intensity"
			)
		#plot(c(xmin,xmax+0.3), c((0.1),(66535)+5000), bg = "black", cex = 0.5, main="data", type="n")
		#points(log(x), log(y[1,]), col="grey", lty=2) 
		
		#now plot the fitted value
		xx<-seq(xmin,xmax, by=(xmax-xmin)/1000)
		type<-match.arg(type);
		#cat("type is :",type,"\n");
		for(i in c(1:dim(y)[1]))
		{
			points(x_shift[i,], y[i,], col=i, cex=0.5
			, type=type
			)
			#yy<-f_PowerBaselinePMT(c(B, G, phi[i], delta),xx);
			#lines(xx,yy,col=i)
		}
		yy<-f_PowerBaselinePMT(c(B, G, phi.ref, delta),xx);
		lines(xx,yy, col=2,lwd=2.5, lty=2)
		if(!is.null(filename)&&(filename!=""))
		{
			dev.off();
		}
	}#end of function
	
	
	
	#'@title function to run fitting for the PMT power model with baseline
	#'@description run 
	#'	fitting for the PMT power model with baseline. The input is a data matrix and
	#'	PMT voltage.
	#'@details It takes in the initial values for parameters, B and G.
	#'	We will estimate phi and delta inside the function. Then fits the data to estimate them.\cr 
	#'assuming \itemize{
	#'		\item{1) data have done background correction} 
	#'		\item{2) data have not been log-transformed.}
	#'		\item{3) data of y must be a datamatrix, x could be a vector.} 
	#'		\item{4) we will estimate the initial values for phi and delta. }
	#'		\item{5) phi is log-transformed in the model. In this way, we force it to be positive.
	#'				otherwise, we might run into trouble of doing log of negative values.}
	#'}
	#' As for the output, this function returns the nls.lm fitted object. It contains the pars for
	#'	B, beta and phi. Phi's are the estimated light influx and they are log-transformed and 
	#'	delta-normalized, which means they phi'/delta. This way we get rid of the individual effect
	#'	of estimated delta, because we assume delta is identical across all data from one PMT scanner.
	##input,
	#y<-log(exp(B)+exp(delta*(log(x)+phi))) <-here we can see the phi<-phi'/delta
	#'@param object data matrix holding the signal intensity data. It is arranged as gene by voltage
	#'	each rows is for one data line to be fitted 
	#'			for a pmt power model with baseline. 
	#'@param	x vector the PMT gain setting vector assuming all the arrays having the identical pmt vector
	#'@param B numeric the initial background level of the signal, assumed to be log-transformed
	#'@param G numeric the maximum possible signal level, 2^16-1 in an 16-bit signal presentation
	# ## ' @param phi vector of initial value of phi. We assume this is log-transformed. We estimate 
	# ## '	one phi for each gene (row). 
	# ## '@param delta numeric initial value of delta, the exponent of the function. It is assumed to be
	# ##'	identical for each array/block.
	#'@param data.aggregated bool indicating whether the data has been aggregated. FALSE by default.
	#'@param residual.mode string indicating whether to do log-transformed residuals or regular ones. 
	#'@param data.name string used to name the output debugging fitting plots
	#'@param minimum.phi smallest possible phi value. Smaller than this is not meaningful.
	#'@return a nls.lm object of the fitting. Note: phi is log-transformed and delta-normalized.
	#'@seealso [f_PowerBaselinePMT()]
	#'@export
	gainAdjust_fit_Pbp_data<-function(object, x, 
								 #data.log=F, 
								 data.aggregated=F,
								 B, G=gainscan:::pmt_saturated, #phi, delta,
								 residual.mode=c("log","regular"), data.name="",
								 minimum.phi=gainscan:::minimum.phi,
								debug=F )
	{
		if(missing(object))
		{
			stop("missing object, please specify....")
		}
		if(missing(x))
		{
			stop("missing input x, please specify....");
		}
		#if(missing(phi))
		#{
		#	stop("missing input phi, please specify....");
		#}
		#if(missing(delta))
		#{
		#	stop("missing input delta, please specify....");
		#}
		#if(missing(delta))
		#{
		#	stop("missing input delta, please specify.....");
		#}
		
		if(class(object)!="matrix")
		{
			stop("object is not of correct format. please specify......")
		}
		if(class(x)!="numeric")
		{
			stop("the input of x is not correct. please check........");
		}
		
		if(class(x)=="numeric")
		{
			 
			if(dim(object)[2]!=length(x))
			{
				stop("the input object matrix and x are not correctly formated, please check........"); 
			}
			#x<-matrix(x, nrow=dim(object)[1], ncol=dim(object)[2], byrow=T);
		}
		if(is.null(data.name)||data.name=="")
		{
			data.name<-as.numeric(format(Sys.time(), "%OS3"))*1000 ###in milliseconds
			data.name<-paste("a",data.name, sep="");
		}
		if(missing(B))
		{
			stop("B is missing, please specify......");
		}
		#data.name<-paste(data.name,".jpg", sep="");
		#do aggregation???
		y<-object;
		if(!data.aggregated)
		{
			y<-gainAdjust.aggregate(object, mode="geometric")
		}
		#now let's estimate initials for phi and delta 
		inits<-gainAdjust_prepareInits_Pbp(y=y, x, aggregated=T, data.log=F)
		delta<-inits$inits.delta;
		phi<-log(inits$inits.phi)
		
		nprint<-0;
		if(debug)
		{
			nprint<-1;
		}
		residual.mode=match.arg(residual.mode);
		
		#now estimate initials.
		
		#now getting data ready to do the fitting......
		#cat("Start doing the fitting ......\n")
		#calling fitting 
		#cat("length(phi):",length(phi),"\n");
		#cat("pars:",length(pars),"\n");
		pars<-c(B,delta, phi);
		#lower<-rep(0,length(pars));
		#lower[c(1,2)]<- -Inf;
		#sink("temp.txt");
		ft<-nls.lm(pars, fn=gainAdjust_fnResPbp, #lower=lower, 
					y=y, x=x, G=G, residual.mode=residual.mode,
					control = nls.lm.control(nprint=nprint, maxit=100)
					);
		flush.console();
		
		#now we need to refit the data individually to refine them
		pars<-ft$par[-c(1,2)];
		
		#fit_phi<-gainAdjust_fit_Pbp_phi(object=y,x=x, B=ft$par[1], delta=ft$par[2], debug=debug)
		#pars<-fit_phi;
		
		pars[pars<minimum.phi]<-minimum.phi;
		ft$par[-c(1,2)]<-pars;
		#sink();
		#do plotting for debugging purpose
		if(debug)
		{
			#now plot the fitting.
			gainAdjust_plot_Pbp_raw(LMfit=ft,y=y, x=x, data.aggregated=T,
						G=G, log="xy",
						filename=NULL);
			gainAdjust_plot_Pbp_raw(LMfit=ft,y=y, x=x, data.aggregated=T,
						G=G, log="xy",
						filename=paste0(data.name, "raw.jpg"));
			#now plot them
			gainAdjust_plot_Pbp(LMfit=ft,y=y, x=x, data.aggregated=T,
						G=G, 
						ref.index=1,  #used to indicate which line is used as the reference, all other ones is shifted towards it.
						log="xy", type="l",
						filename=NULL);
			gainAdjust_plot_Pbp(LMfit=ft,y=y, x=x, data.aggregated=T,
						G=G, 
						ref.index=1,  #used to indicate which line is used as the reference, all other ones is shifted towards it.
						log="xy", type="l",
						filename=paste0(data.name, "shift.jpg"));
		}
		
		return(ft);
	
	}#end of the function
	
	#'@title function to run fitting for the PMT power model with baseline
	#'@description run 
	#'	fitting for the PMT power model with baseline to esimtate the individual phi's based on a fixed model (fixed B, delta and G). The input is a data matrix and
	#'	PMT voltage.
	#'@details It takes in the parameters, B, phi and delta, which are estimated by fitting previously.
	#'	Then fits the data to estimate phi's for each gene data series. This usually is done for expression data based on the control fitted model.\cr 
	#'assuming \itemize{
	#'		\item{1) data have done background correction} 
	#'		\item{2) data have not log'ed, but aggregated or don't care whether it is aggregated.
	#'					it the caller's responsibility to make it is aggregated or not.}
	#'		\item{3) data of y must be a datamatrix, x could be a vector.} 
	#'		\item{4) we have the initials estimated inside the function. }
	#'		\item{5) phi is log-transformed. In this way, we force it to be positive.
	#'				otherwise, we might run into trouble of doing log of negative values.}
	#'}
	#	#' As for the output, this function returns estimated phi. 
	#'	 Phi's are the estimated light influx and they are log-transformed and 
	#'	delta-normalized, which means they phi'/delta. This way we get rid of the individual effect
	#'	of estimated delta, because we assume delta is identical across all data from one PMT scanner.
	##input,
	#y<-log(exp(B)+exp(delta*(log(x)+phi))) <-here we can see the phi<-phi'/delta
	#'@param object data matrix holding the signal intensity data. It is arranged as gene by voltage
	#'	each rows is for one data line to be fitted 
	#'			for a pmt power model with baseline. 
	#'@param	x vector the PMT gain setting vector assuming all the arrays having the identical pmt vector
	#'@param B numeric the initial background level of the signal. assumed to be log-transformed
	#'@param G numeric the maximum possible signal level, 2^16-1 in an 16-bit signal presentation
	# #' ## @param phi vector of initial value of phi. We assume this is log-transformed. We estimate 
	# #'	one phi for each gene (row). 
	#'@param delta numeric initial value of delta, the exponent of the function. It is assumed to be
	#'	identical for each array/block.
	#'@param residual.mode string indicating whether to do log-transformed residual or regular. 
	# # '@param minimum.phi a numeric number indicating the lowest possible phi. Smaller than this is not meaningful.
	#'@return a nls.lm object of the fitting. Note: phi is log-transformed and delta nomalized.
	#'@seealso [f_PowerBaselinePMT()]
	#'@export
	gainAdjust_fit_Pbp_phi<-function(object, x, 
								#ylong=F, 
								# data.aggregated=F,
								 B, G=gainscan:::pmt_saturated, delta,
								 residual.mode=c("log","regular"),
								 #minimum.phi=minimum.phi,
								debug=F )
	{
		if(missing(object))
		{
			stop("missing object, please specify....")
		}
		if(missing(x))
		{
			stop("missing input x, please specify....");
		}
		if(missing(delta))
		{
			stop("missing input delta, please specify....");
		}
		
		if(class(object)!="matrix")
		{
			stop("object is not of correct format. please specify......")
		}
		if(class(x)!="numeric")
		{
			stop("the input of x is not correct. please check........");
		}
		
		if(class(x)=="numeric")
		{
			 
			if(dim(object)[2]!=length(x))
			{
				stop("the input object matrix and x are not correctly formated, please check........"); 
			}
			#x<-matrix(x, nrow=dim(object)[1], ncol=dim(object)[2], byrow=T);
		}
		
		nprint<-0;
		if(debug)
		{
			nprint<-1;
		}
		residual.mode=match.arg(residual.mode);
		#now getting data ready to do the fitting......
		#cat("Start doing the fitting ......\n")
		#calling fitting 
		#cat("length(phi):",length(phi),"\n");
		#cat("pars:",length(pars),"\n");
		
		
		inits<-gainAdjust_prepareInits_Pbp(y=object, x, aggregated=T, data.log=F);
		#inits$inits.phi[inits$inits.phi>0.01]<-median(inits$inits.phi)
		pars<-log(inits$inits.phi);
		phi<-pars; #get ready for output.
		#pars<-c(phi);
		#lower<-rep(0,length(pars));
		#lower[c(1,2)]<- -Inf;
		for(i in 1:dim(object)[1]){
			#cat("doing.....i:",i,";par[i]:",pars[i],"\n");
			ft<-nls.lm(pars[i], fn=gainAdjust_fnResPbp_phi, #lower=lower, 
						y=as.numeric(object[i,]), x=x, B=B, G=G, delta=delta, residual.mode=residual.mode,
						control = nls.lm.control(nprint=nprint, maxit=100)
						);
			phi[i]<-ft$par;
		}
		phi[phi<minimum.phi]<-minimum.phi;
		return(phi);
	
	}#end of the function

	
#'@title function to run fitting for the PMT power model with baseline
	#'@description this is the function to 
	#'	fit one array for the PMT power model with baseline. The input is an array (EListRaw 
	#'	object {\code{\link{EListRaw-class}}}.
	#'@details It takes in the initial values for parameters G. The initial values 
	#'	for B, phi and delta are estimated internally.
	#' Then fits the array data to estimate them. 
	# # '	blocks and also aggregate them.  \cr 
	#'As for the output, this function returns list. It contains the EListRaw object for 
	#'	phis and B and delta estimated as well. Phi's are the estimated light influx and they are log-transformed and 
	#'	delta-normalized, which means they phi'/delta. This way we get rid of the individual effect
	#'	of estimated delta, because we assume delta is identical across all data from one PMT scanner.
	##input,
	#y<-log(exp(B)+exp(delta*(log(x)+phi))) <-here we can see the phi<-phi'/delta
	#'@param object EListRaw object array data. Inside this object, the data matrix
	#'	tables should arranged as gene (row) by voltage (column). The columns are 
	#'	orders as x. 
	#'@param	x vector the PMT gain setting vector assuming all the arrays having the identical pmt vector
	#'
	#'@param G numeric the maximum possible signal level, 2^16-1 in an 16-bit signal presentation
	#'@param block.size number of blocks to used for one fitting.
	#'@param fit.mode text to indicate whether to control or all data (control+Expression) data for fitting.
	#'@param data.aggregated bool to indicate whether the data has been aggregated.
	#'@param residual.mode string indicating whether to do log-transformed residual or regular. 
	#'@param debug bool to indicate whether do debugging mode, by plotting and writing extra information
	#'@return a data list (EListRaw object) contains following the fitted phi's for E and C. 
	#'	Note: the phi's returned are log-transformed,
	#'	and also was delta-normalized.  
	#'@seealso [f_PowerBaselinePMT()] [gainAdjust_fit_Pbp_data()]
	#'@export
	gainAdjust_fit_Pbp_array<-function(object, x, 
								#ylong=F, 
								data.aggregated=F,
								 B=gainscan::pmt_baseline, 
								G=gainscan::pmt_saturated, #phi, delta,
								block.size=4, fit.mode=c("control","all"),
								 residual.mode=c("log","regular"), array.name="",
								debug=F )
	{
		if(missing(object))
		{
			stop("missing object, please specify....")
		}
		if(missing(x))
		{
			stop("missing input x, please specify....");
		}
		
		if(class(object)!="EListRaw"&&class(object)!="list")
		{
			stop("object is not of correct format. please specify......")
		}
		if(class(x)!="numeric")
		{
			stop("the input of x is not correct. please check........");
			if(dim(object$E)[2]!=length(x)||dim(object$C)[2]!=length(x))
			{
				stop("the input object matrix and x are not correctly formated, please check........"); 
			}
			#x<-matrix(x, nrow=dim(object)[1], ncol=dim(object)[2], byrow=T);
		}
		ydata<-object;
		if(!data.aggregated)
		{
			ydata<-gainAdjust.aggregate(object, mode="geometric")
		}
		if(is.null(array.name)||array.name=="")
		{

			array.name<-as.numeric(format(Sys.time(), "%OS3"))*1000 ###in milliseconds
			array.name<-paste("a",array.name, sep="");
		}
		#paste0(data.name, "raw.jpg")
		
		fit.mode<-match.arg(fit.mode);
		residual.mode<-match.arg(residual.mode);
		blkNum<-length(unique(ydata$gene[,"Block"]))
		blocks<-unique(ydata$gene[,"Block"]);
		#block.size<-match.arg(block.size, choices=c(1:max(blkNum)))
		
		#now getting data ready to do the fitting......
		#get output ready
		phi<-ydata; #a holder for
		phi$C<-matrix(ydata$C[,1],ncol=1);
		phi$E<-matrix(ydata$E[,1],ncol=1);
		phi$targets<-ydata$targets[1,];
		#fits<-vector(mode="list", length=ceiling(blkNum/block.size));
		Bs <-rep(0, 	ceiling(blkNum/block.size));
		deltas<-rep(0,ceiling(blkNum/block.size));
		for(j in 1:ceiling(blkNum/block.size)){
			cat("\tFitting data by blocks: ",j, "/",blkNum/block.size, "...\n")
			#for the array, we need to get the protein data by block, 
			bg<-(j-1)*block.size+1
			ed<-(j-1)*block.size+block.size
			#cat(paste("bg:",bg, "; ed:", ed,"\n", sep=""));
			if(ed>blkNum)
			{
				ed<-blkNum
			}
			f.name =paste0(array.name, "_block_", j);#cat("file anme is", f.name,"\n");
			idx.prC<-which(is.element(ydata$cgene[,"Block"],blocks[c(bg:ed)]))
			ydata.blocked<-ydata$C[idx.prC,]
			idx.prE<-which(is.element(ydata$gene[,"Block"],blocks[c(bg:ed)]))
			if(fit.mode=="all")
			{
				ydata.blocked<-rbind(ydata$E[idx.prE,], ydata.blocked)
			}
			#now data are ready, we need to do aggregration/inits
			#inits<-gainAdjust_prepareInits_Pbp(y=ydata, x, aggregated=F, data.log=F)
			#B<-log(1.1)
			
			#sink("temp.txt");
			#system.time(
			ftf<-gainAdjust_fit_Pbp_data(object=ydata.blocked, x=x, 
											#ylong=F, 
											 data.aggregated=T,
											 B=B, G=G, #phi=log(inits$inits.phi), delta=inits$inits.delta,
											 residual.mode="log", data.name=f.name,
											debug=debug)
					#		)
			#sink();
			#now we got the fitting, so we need to write out the output.
			Bs[j]<-ftf$par[1]
			deltas[j]<-ftf$par[2];
			B<-ftf$par[1];
			delta<-ftf$par[2];
			ftf$par<-ftf$par[-c(1,2)]
			
			#write output of expression targets.
			if(fit.mode=="all")
			{
				phi$E[idx.prE,1]<-ftf$par[c(1:length(idx.prE))]
				#if(j==1)
				#{
				#	phi$E<-ftf$par[c(1:length(idx.prE))]
				#}else{
				#	phi$E<-rbind(phi$E,ftf$par[c(1:length(idx.prE))]); #we did not changed the order.
				#}
				ftf$par<-ftf$par[-c(1:length(index.prE))];
			} else #need to fit using the parameters based on the controls.
			{
				#cat("\nDoing fit for E proteins.......\n");
				ftf_single<-gainAdjust_fit_Pbp_phi(object=ydata$E[idx.prE,], x=x, B=B,
								G=G, delta=delta, residual.mode=residual.mode,debug=debug);
				phi$E[idx.prE,1]<-ftf_single;
				#if(j==1)
				#{
				#	phi$E[idex.prE,1]<-ftf_single;
				#}else{
				#	phi$E<-rbind(phi$E, matrix(ftf_single,ncol=1));
				#}
			}
			#write output for control targets
			phi$C[idx.prC,1]<-ftf$par;
			#if(j==1)
			#{
			#	phi$C[idx.prC,1]<-ftf$par
			#}else
			#{
			#	phi$C<-rbind(phi$C,matrix(ftf$par,ncol=1));
			#}	
			cat("\n");
			flush.console();
		}#for loop for different blocks of the object
		return(list("elist"=phi, "block.size"=block.size, "B"=matrix(Bs, ncol=1),"delta"=matrix(deltas, ncol=1)));
	}#End of the function fit array
	
	#'@title function to fitting array data the PMT power model with baseline 
	#'@description this is the function to 
	#'	fit array data with the PMT power model with baseline. The input is protoArray (EListRaw 
	#'	object) {\code{\link{EListRaw-class}}}.
	#'@details It takes in the initial values for parameters G. The initial values 
	#'	for B, phi and delta are estimated internally.
	#' Then fits the array data to estimate them. 
	#'  We read the EListRaw$targets table to look for different arrays (filed of "Array")
	#'	. Then group the data from the same array together. Within each array data group, we assume
	#'	the data are arranged to have the same order of PMT voltages (the vector of x, see below).
	#'As for the output, this function returns the list. It contains the pars for
	#'	B, beta and EListRaw. The EListRaw object contains Phi's,
	#	which are the estimated light influx and they are log-transformed and 
	#'	delta-normalized, which means they phi'/delta. This way we get rid of the individual effect
	#'	of estimated delta, because we assume delta is identical across all data from one PMT scanner.
	##input
	#'@param object EListRaw object holding data for all the arrays. Inside this object, the data matrix
	#'	tables should arranged as gene (row) by voltage (column). The columns are 
	#'	orders as x. It contains the data for different arrays and each arrays are read under 
	#'	different PMT voltages. The design information is included in the table "targets".
	#'@param	x vector the PMT gain setting vector assuming all the arrays having the identical pmt vector
	#'
	#'@param G numeric the maximum possible signal level, 2^16-1 in an 16-bit signal presentation
	#'@param block.size number of blocks to used for one array fitting.
	#'@param fit.mode text to indicate whether to control or all data (control+Expression) data for fitting.
	#'
	#'@param residual.mode string indicating whether to do log-transformed residual or regular. 
	#'@param debug bool to indicate whether do debugging mode, by plotting and writing extra information
	#'@return a data list with the following elements:
	#' elist:	(EListRaw object) contains the fitted phis for E and C. 
	#'	Note: the phi's returned are log-transformed. To get the phi as the light influx, we need
	#'	to transform them back by (exp(phi))^delta. 
	#'B: fitted B. The number of elements depends on the block size. It arranged as block by array. 
	#'	each array is one column and number of element depends on the block.size.
	#'delta: fitted delta. The number of elements depends on the block.size. It arranged as block by array. 
	#'	each array is one column and number of element depends on the block.size.
	#'block.size: number of block.size to fit for each array data. 
	#'@seealso [f_PowerBaselinePMT()] [gainAdjust_fit_Pbp_data()] [EListRaw-class]
	#'@export
	gainAdjust_fit_Pbp<-function(object, x, 
								#ylong=F, aggregated=F,
								 B=gainscan::pmt_baseline, 
								G=gainscan::pmt_saturated, #phi, delta,
								block.size=4, fit.mode=c("control","all"),
								 residual.mode=c("log","regular"),
								debug=F )
	{
		if(missing(object))
		{
			stop("missing object, please specify....")
		}
		if(missing(x))
		{
			stop("missing input x, please specify....");
		}
		
		if(class(object)!="EListRaw"&&class(object)!="list")
		{
			stop("object is not of correct format. please specify......")
		}
		if(class(x)!="numeric")
		{
			stop("the input of x is not correct. please check........");
			if(dim(object$E)[2]!=length(x)||dim(object$C)[2]!=length(x))
			{
				stop("the input object matrix and x are not correctly formated, please check........"); 
			}
			#x<-matrix(x, nrow=dim(object)[1], ncol=dim(object)[2], byrow=T);
		}
		#aggregate 
		ydata<-gainAdjust.aggregate(object, mode="geometric")
		
		#parse the data into single array and then call to do fitting.
		arrays<-unique(ydata$targets[,"Array"])
		phi<-c();
		Bs<-c();
		deltas<-c();
		block.sizes<-c();
		cat("Starting fitting the array data with the PMT power model.........\n");
		for(i in 1:length(arrays))
		{
			cat("Fitting array ", i, "/", length(arrays), "....\n");
			array.index<-grep(arrays[i],ydata$targets[,"Array"])
			ydata.array<-ydata;
			ydata.array$E<-ydata$E[,array.index]
			ydata.array$C<-ydata$C[,array.index]
			ydata.array$targets<-ydata$targets[array.index[1],];
			
			#now feed it to the fitting function
			phi_single<-gainAdjust_fit_Pbp_array(object=ydata.array, x=x, B=B, 
								G=G, data.aggregated=T,
								block.size=block.size, fit.mode=fit.mode,
								 residual.mode=residual.mode, array.name=arrays[i],
								debug=debug );
			#now put things together
			if(i==1)
			{
				phi<-phi_single$elist ;
				Bs<-phi_single$B
				deltas<-phi_single$delta
				
			} else {
				phi$E<-cbind(phi$E,phi_single$elist$E);
				phi$C<-cbind(phi$C,phi_single$elist$C);
				phi$targets<-rbind(phi$targets, phi_single$elist$targets);
				
				Bs<-cbind(Bs, phi_single$B)
				deltas<-cbind(deltas, phi_single$delta)
			}
			block.sizes<-c(block.sizes, phi_single$block.size)
		}##end of for loop
		colnames(phi$E)<-paste0("a",phi$targets[,"Array"])
		colnames(phi$C)<-paste0("a",phi$targets[,"Array"])
		cat("Done..........\n");
		#phi$E<-phi$E*mean(deltas);
		#phi$C<-phi$C*mean(deltas);
		return(list("elist"=phi, "block.size"=block.sizes, "B"=Bs, "delta"=deltas));
	}#End of the function fit array
	
	#'@title scale up the estimated light influx by delta
	#'@description scale up the estimated light influx by the estimated delta parameter.
	#' the esimtated light influx, phi, is log-transformed and delta-normalized. This way
	#' the effect of each individual estimation block is removed. So to get back to its original
	#' scale we need to scale/multiple delta to it. We also need to use the mean delta across all
	#'	estimations using the same PMT to do the scaling, because we assume it is the same as long
	#' as the same PMT is used. Without scaling, it doesn't change the comparison of variance, etc,
	#'	but does affect the scale/magnitude.
	#'@param object EListRaw object that containing the log-transformed fitted phis
	#'@param delta numeric that is the mean of delta for all the estimates generated with one PMT
	#'@return an EListRaw object scaled/multipled by mean delta.
	#'@export
	
	scaleByDelta<-function(object, delta)
	{
	
		if(missing(object)||class(object)!="EListRaw")
		{
			stop("please specify the input object as EListRaw");
		}
		if(missing(delta))
		{
			stop("Please specify the input delta as a numeric");
		}
		object$E<-object$E*delta
		object$C<-object$C*delta;
		#we are done. don't touch anything else.
		return(object)
	}
	
	#'@title unscale the phi objects by delta
	#'@description this is the reverse functionof [scaleByDelta()] by removing the scale 
	#'	of delta.
	#'@export
	unScaleByDelta<-function(object, delta)
	{
	
		if(missing(object)||class(object)!="EListRaw")
		{
			stop("please specify the input object as EListRaw");
		}
		if(missing(delta))
		{
			stop("Please specify the input delta as a numeric");
		}
		object$E<-object$E/delta
		object$C<-object$C/delta;
		#we are done. don't touch anything else.
		return(object)
	}
	#'@title acquire the feature spot flourescence
	#'@description calculate the feature spot flourescence intensity based on the 
	#' power function model, as if we are scan the array to acquire the data using 
	#' a scanner.
	#'@details The following equation is used to convert/amplify the input signal to
	#'	the flourescence intensity in the protein array\cr
	#'	log(Y_i)=log(B+phi_i*(v)^delta)+epsilon_ij \cr
	#'	epsilon_ij ~ N(0, sigma^2)
	#'@param object list of objects return by [gainAdjust_fit_Pbp()] containing \cr
	#'	\itemize{
	#'		\item(1)elist this is the phi's (input signals) for all the spots on array
	#'		
	#'		\item(2)B: fitted B. The number of elements depends on the block size. 
	#'				It arranged as block by array. 
	#'				each array is one column and number of element depends on the block.size.
	#'		\item(3)delta: fitted delta. The number of elements depends on the block.size. 
	#'				It arranged as block by array. 
	#'				each array is one column and number of element depends on the block.size.
	#'		\item(4)block.size: number of block.size to fit for each array data.
	#'	}
	#'@param numeric the voltage to acquire the data
	#'@param normalized TRUE means all parameters are fixed for all fitting blocks. By 
	#'	default (False) we will used fitting block specific parameters to get best fitted lines
	#'@return an EListRaw objects with the 
	#'@export
	#<----------+++now we are doing individual parameter (no normalizatoin)
	# # will try the normalization method later.......
	acquire<-function(object, voltage, PMT.Saturated=gainscan:::pmt_saturated, normalized=FALSE)
	{
		if(missing(object))
		{
			stop("missing input object, please specify!!!")
		}
		if(missing(voltage))
		{
			stop("missing input voltage, please specify")
		}
		if(class(object)!="list"||length(object)!=4)
		{
			stop("the input object is not correctly constructed, please check....");
		}
		if(length(voltage)>1)
		{
			voltage<-voltage[1];
			warning("More than one elements were specified for voltage. only the first one will be used. ")
		}
		
		#make the parameters ready
		array.num<-dim(object$elist$targets)[1];
		if(array.num!=length(object$block.size))
		{
			stop("the block.size and array don't have the identical size, please double check!!!")
		}
		object.ret<-object$elist;
		
		#block.size<-object$block.size;
		#max.block.size<-sum(block.size)
		for(j in 1:array.num)
		{
			
			block.size<-object$block.size[j];
			max.block.num<-max(object$elist$gene$Block)
			for(i in 1:ceiling(max.block.num/block.size))
			{
				###now need to figure out the smaple blocks
				block.size.start<-block.size*(i-1)+1
				block.size.end<-block.size.start+block.size-1;
				if(block.size.end > max.block.num)
				{
					block.size.end<-max.block.num;
				}
				index<-which(object$elist$genes$Block>= block.size.start 
							& object$elist$genes$Block<= block.size.end);
				#now get the parameter
				B<-object$B[i,j];
				delta<-object$delta[i,j];
				
				#now depending on the mode (normalized vs un-normalized)
				#we do acquisition differently
				if(normalized)
				{
					B<-mean(object$B)
					delta<-mean(object$delta)
					object.ret$E[index,j]<-
						f_PowerBaselinePMT_phi(B=B,G=PMT.Saturated , phi=object$elist$E[index,j], delta=delta, x=voltage, saturation=F);
					#for control, at this moment we assume control and target has the same set of parameters
					index<-which(object$elist$cgenes$Block>= block.size.start 
							& object$elist$cgenes$Block<= block.size.end)
					object.ret$C[index,j]<-
						f_PowerBaselinePMT_phi(B=B,G=PMT.Saturated , phi=object$elist$C[index,j], delta=delta, x=voltage, saturation=F);
					cat("doing normalized mode\n")
					next;
				}	
				cat("doing unnormalized mode\n")
				object.ret$E[index,j]<-
					f_PowerBaselinePMT_phi(B=B,G=PMT.Saturated , phi=object$elist$E[index,j], delta=delta, x=voltage, saturation=F);
				#for control, at this moment we assume control and target has the same set of parameters
				index<-which(object$elist$cgenes$Block>= block.size.start 
							& object$elist$cgenes$Block<= block.size.end)
				object.ret$C[index,j]<-
					f_PowerBaselinePMT_phi(B=B,G=PMT.Saturated , phi=object$elist$C[index,j], delta=delta, x=voltage, saturation=F);
			}	
		}
		return(object.ret);
	}#end of function
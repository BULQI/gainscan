#####Description for this file runStats.R#####
# this file is used to define the functions to run statistical analysis 
# in order to identify the proteins differentially binding to the testing Abs
############### 

# the following section is used to act like a holder for
# importing libraryies
# package limma is necessary for EList
# # ' @title an empty file to import necessary libraries 
#' @import limma 
## ' @name _empty
NULL

#first define the top most function to drive the analysis
#the overall work flow
#'@title S3 function to identify features with significant interaction
#'@description \code{identifyFeatures} to identify features/proteins that interact 
#'	siginificantly with the testing Abs
#'@details This function takes in an elist object holding the correctly
#'	 normalized protoarray data, and then run statistical analysis 
#' 	 using an linear regression model with interaction terms. The results 
#'	 are also corrected for multiple analysis using FDR
#'	(see details about the model and analysis {\code{\link{gainscan}}})
#'	Please see the example below for function usage.
#'
#'@param dataFilePath string path to the file containing the text formatted input data.
#'				This is a tab-delimited file exported from the software reading
#'               the arrays. It has both gene data and control data in one.
#'@param TargetFile string path to a text file holding the meta data for the array. 
#'				Please check the template
#'@param start.data numeric the line number indicating where the data section starts
#'@param nrows.data numeric number of rows for the data section
#'@param start.control numeric the line number indicating where the control data section starts
#'@param nrow.control numeric the number of rows for the control section
#'@param aggregation string indicating which type of duplicate aggregation should be performed. 
#			If "min" is chosen, the value for the corresponding feature will be the minimum of both. 
#			If "arithMean" is chosen, the arithmetic mean will be computed. 
#			If "genMean" is chosen, the geometric mean will be computed.
#			The default is "min" (optional).
#'@param as.is boolean TRUE by default,for reading the text file
#'@param header boolean TRUE by default, for reading the text file
#'@param sep char '	' (tab) by default, for reading the text file
#'@param na.strings string "" by default, for reading the text file
#'@param quote char "" (no quote) by default,  for reading the text file
#'@return an ELISTlist object containing both the data and control
#@seealso
#'@examples
#'	datapath<-system.file("extdata", package="ARPPA")
#'	targets <- list.files(system.file("extdata", package="ARPPA"),
#'		 pattern = "targets_text_Batch1", full.names=TRUE) 

#'	elist2<-importTextData(dataFilePath=datapath, targetFile=targets, start.data=51, nrows.data=18803-53,
#'				start.control=18803, nrows.control=23286,aggregation="geoMean",
#'				as.is=TRUE, header=TRUE,sep="\t",na.strings="", quote="")
#'
#'@export
identifyFeatures<-function(elist, expressionCutOff=50 )
		{
		}
		
#calculate the fisher signal to noise ratio
#      s=(u1+u2)^2/(var1+var2)
#  (definition based on sboner et al work).
#'@title fisher's signal-to-noise ratio
#'@description calculate the fisher's signal-to-noise ratio between two groups
#'@details This works on two EListRaw objects. Each object contains one group 
#'	of arrays.  
#'@param object1 EListRaw object holding the data for one group of arrays
#'@param object2 EListRaw object holding the data for one group of arrays
#'@param mode indicating different way of calculating the SNR
#'	var: snr=(u1-u2)^2/(s1^2+s2^2)
#'	std: snr=abs(u1-u2)/(s1+s2)
# #'	sqrt: snr=sqrt(snr_var) (see the first one named "var")	
#'@return a list containing two fields. One fields of fisher signal to noise ratio for 
#'	feature proteins and another for control proteins.
#'@export
#
signalToNoiseRatio<-function(object1, object2, mode=c("var","std")#,"sqrt")
							)
{
	if(class(object1)!="EListRaw"||class(object2)!="EListRaw")
	{
		stop("input is not in a correct format, please check!!!");
	}
	if(dim(object1$E)[2]<=1||dim(object1$C)[2]<=1||dim(object2$E)[2]<=1||dim(object2$C)[2]<=1)
	{
		stop("input doesn't has replicates, please check");
	}
	mode<-match.arg(mode);
	#doing E fields first
	us<-apply(object1$E,1,mean)-apply(object2$E,1,mean);
	vars<-apply(object1$E,1,var)+apply(object2$E,1,var);
	vars[vars==0]<-0.1
	E<-NULL;
	C<-NULL;
	switch(mode,
		var={E<-us*us/vars;},
		std={E<-abs(us)/(sqrt(apply(object1$E,1,var))+sqrt(apply(object2$E,1,var)));
			},
		#sqrt={E<-sqrt(us*us/vars);}
	)
	
	us<-apply(object1$C,1,mean)-apply(object2$C,1,mean);
	vars<-apply(object1$C,1,var)+apply(object2$C,1,var);
	vars[vars==0]<-0.1
	switch(mode,
		var={C<-us*us/vars;},
		std={C<-abs(us)/(sqrt(apply(object1$C,1,var))+sqrt(apply(object2$C,1,var)));
			},
		#sqrt={C<-sqrt(us*us/vars);}
	)
	#C<-us*us/vars;
	
	s<-list(E=E,C=C);
}
#calling to do t.tests in batch. and return 
#dtf either a data.frame or matrix
#assuming the groups are by columns. We will run t tests by row and separate columns 
#into groups by ind.x and ind.y (no overlapping)
#'@title running t test in batch 
#'@description calling t.test to run then in batch to process
#'	data from data.frame 
#'@export
ttests_batch<-function(dtf1,dtf2, alternative = c("two.sided", "less", "greater"))						
{
	if((class(dtf1)!="data.frame"&&class(dtf1)!="matrix")||(class(dtf2)!="data.frame"&&class(dtf2)!="matrix"))
	{
		stop("input is not in a correct format, please check!!!");
	}
	#if(!all(!is.element(ind.x, ind.y)))
	#{
	#	stop("input indices have overlapping, please check");
	#}
	num<-min(dim(dtf1)[1],dim(dtf2)[1])
	sts<-rep(0, dim(dtf1)[1]);
	ps<-rep(0, dim(dtf1)[1]);
	for(i in 1:num)
	{
		if(var(dtf1[i,])==0)
		{
			dtf1[i,]<-dtf1[i,]+rnorm(length(dtf1[i,]),0,1)
		}
		if(var(dtf2[i,])==0)
		{
			dtf2[i,]<-dtf2[i,]+rnorm(length(dtf2[i,]),0,1)
		}
		tts<-t.test(dtf1[i,],dtf2[i,],alternative,var.equal=T)#t.test(dtf1[i,],dtf2[i,],alternative)
		sts[i]<-tts$statistic
		ps[i]<-tts$p.value
	}
	
	list(statistic=sts, p.value=ps);
}

#calling to do t.tests in batch. and return 
#dtf either a data.frame or matrix
#assuming the groups are by columns. We will run t tests by row and separate columns 
#into groups by ind.x and ind.y (no overlapping)
#'@title running t test on ELIST_raw 
#'@description calling t.test to process
#'	data from EListRaw objects 
#'@export
ttests<-function(obj1, obj2, alternative = c("two.sided", "less", "greater"))						
{
	if((class(obj1)!="EListRaw")||(class(obj2)!="EListRaw"))
	{
		stop("input is not in a correct format, please check!!!");
	}
	#if(!all(!is.element(ind.x, ind.y)))
	#{
	#	stop("input indices have overlapping, please check");
	#}
	#num<-min(dim(obj1)[1],dim(dtf2)[1])
	#sts<-rep(0, dim(dtf1)[1]);
	#ps<-rep(0, dim(dtf1)[1]);
	#for(i in 1:num)
	#{
	#	tts<-t.test(dtf1[i,],dtf2[i,])
	#	sts[i]<-tts$statistic
	#	ps[i]<-tts$p.value
	#}
	tts.E<-ttests_batch(obj1$E,obj2$E,alternative);
	tts.C<-ttests_batch(obj1$C,obj2$C, alternative);
	list(E=tts.E, C=tts.C);
}
#'@title plot SNR vs expression level
#'@description plot SNR vs the expression level
#'@details TBD
#'@param SNR.mode character for mode to calculate SNRs.
#'@param add boolean to indicate wheter to add to the current graph
#'@param col integer to add col to plotting.
#'@param pch integer to plot the points
#'@param log characters to indicate whether to plot (add=T) 
#'		the figures in log scale for axes.
#'@param plot.mode characters to indicate whether to plot mean expression or difference
#'		between 2 groups vs SNR
#'		mean, mean expression between two groups
#'		difference, the abolute level of group expression difference.
#'@export
plotSNR<-function(object1, object2, SNR.mode=c("var","std"),
					plot.mode=c("mean","difference"),
					add=F,col=1, pch=1,log="y"
							)
{
	#cat("begining")
	if(class(object1)!="EListRaw"||class(object2)!="EListRaw")
	{
		stop("input is not in a correct format, please check!!!");
	}
	if(dim(object1$E)[2]<=1||dim(object1$C)[2]<=1||dim(object2$E)[2]<=1||dim(object2$C)[2]<=1)
	{
		stop("input doesn't has replicates, please check");
	}
	#cat("start doing it");
	SNR.mode<-match.arg(SNR.mode);
	plot.mode<-match.arg(plot.mode);
	snr<-signalToNoiseRatio(object1, object2, SNR.mode)
	#now calculate the mean of expression for E only
	meanExp<-NULL;
	if(plot.mode=="mean"){
		meanExp<-(apply(object1$E,1,mean)+apply(object2$E,1,mean))/2
	} else if(plot.mode=="difference"){
		meanExp<-abs(apply(object1$E,1,mean)-apply(object2$E,1,mean));
	}
	#print("feng")
	#cat("before")
	if(!add){
		plot(meanExp,snr$E, xlab="mean expression", ylab="SNR", main="expression vs SNR",log=log, col=col,pch=pch)
	}else {
		points(meanExp,snr$E, xlab="mean expression", ylab="SNR", main="expression vs SNR",
				col=col,pch=pch)
		}
	cat("done!!\n")
#	return (1);
}

#regression by matrix to estimate linear parameters
#"manually" do this instead of lm to save (?) some computer time
#hoping the matrix is not too huge 
#'@export

regressionByMatrix<-function(y, x)
{
	if(class(y)!="numeric"||class(x)!="numeric")
	{
		stop("please specify the input (x and y) as matrix")
	}
	if(length(y)!=length(x))
	{
		stop("x and y are of not equal length, please check")

	}
	
	#build the matrix
	X<-matrix(1, nrow=length(x))
	X<-cbind(X, x)
	Y<-matrix(y, nrow=length(y))
	
	bt<-solve(crossprod(X))%*%crossprod(X,Y)
}
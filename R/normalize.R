###this is the R code to do the protoArray normalization (Sboner et al)
## with the capacity of doing 5PL fitting
##started 1/2/2018
##---Feng @BU#####

#'@include data.R 
#used to compile code from data.R file, because we need the function
#eg. matrix2dframe function

#this is the 
#normalizeArrays<-function()
#{
#
#}

#this is the code to run normalization using RLM only (no 5PL yet)
#'@title normalize array using the Robust Linear Model
#'@description implementation of the Robust Linear Model following Sboner, et al work.
#'@details Only RLM implemented so far.
#'@return the normalized the array data in EListRaw format. 
#	
#'@export

normalizeArrays.RLM<-function(data,controls="HumanIg",
				method=c("RLM"), #only implemented RLM for now
				log=TRUE, log.base=exp(1), coding=c("Deviation","Dummy","Simple")
			)
{
	#check the input 
	if(missing(data))
	{
		stop("***Error: input data missing, please specify!!");
	}
	if(class(data)!="EListRaw")
	{
		stop("***Error: unknown input data format, please use EListRaw data"); 
	}
	if(missing(coding))
	{
		stop("**** please specify the coding for RLM fitting"); 
	}
	method<-match.arg(method)
	coding<-match.arg(coding)
	if(coding!="Deviation")
	{
		msg<-"WARNING: the coding is chosen to be Dummy or Simple, which are not carefully debugged. Please select to do Deviation instead!!!";
		warning(msg)
		stop(msg);
	}
	#now everything seems OK
	#get the data set ready
	if(dim(data$C)[2]!=dim(data$E)[2])
	{
		stop("Error: the control and feature data contain unequal number of arrays, please check!!")
	}
	if(dim(data$E)[2]<2||dim(data$C)[2]<2)
	{
		stop("Error: the control/feature data contain single array, no normalization is necessary!!")
	}
	#cat("1\n");
	ndata<-data; #ndata will be returned
	if(log){
		ndata$E<-log(data$E, base=log.base);
		ndata$C<-log(data$C, base=log.base);
	}
		#cat("2\n");
	####now get the data ready for linear fitting
	dtf.C<-cmatrix2dataframe(ndata$C, ndata$cgenes, control=controls, log=FALSE);
		#cat("3\n");
	#do fitting to estimate the parameters, array and block, assuming a linear model no interaction
	fit.RLM.C<-fit.RLM(dtf.C, coding=coding);
		#cat("4\n");
	#now need to parse the estimated parameters
	factor.levels<-c(dim(data$C)[2],length(unique(data$cgenes$Block)))
	#factor.first<-c(colnames(data$C)[1],unique(data$cgenes$Block)[1])
		#cat("5\n");
	#factor.first<-c(0,0)
	#arrayNames<-rownames(data$C);
	#if(is.null(arrayNames))
	#{
	#	factor.first[1]<-"1";
	#} else {
	#	facor.first[1]<-arrayNames[1];
	#}
	#factor.first[2]<-data$cgenes$Block[1];
	fct.list<-readFactors(fit.RLM.C, factors=c("array","block"),factor.levels);#,factor.first);
	
	#now can do the normalization by subtracting the factor effects
#	if(coding=="Dummy"){
#		ndata<-sub.effects.Dummy(ndata,factor.effects=fct.list)
#	} else {
#		ndata<-sub.effects.Simple(ndata,factor.effects=fct.list)
#	}	
	ndata<-sub.effects(ndata,factor.effects=fct.list, coding="Deviation");
	
	#----first doing normalization of E field (features)
	#+arrayNames<-colnames(data$E);
	#+blockNames<-(data$genes$Block);
		#cat("6\n");
	#using the first array and block as the reference, and then normalize the rest array and block
	#array effects
	#+fcts<-matrix(rep(fct.list[["array"]][arrayNames],dim(data$E)[1]),
	#+			nrow=dim(data$E)[1], ncol=dim(data$E)[2],byrow=T);
	#+ndata$E<-ndata$E-fcts;
	
	#block effects
	#+fcts<-matrix(rep(fct.list[["block"]][blockNames],dim(data$E)[2]),
	#+			nrow=dim(data$E)[1], ncol=dim(data$E)[2],byrow=F);
	#+ndata$E<-ndata$E-fcts;
	
	#-----now doing the control too
	#+arrayNames<-colnames(data$C);
	#+blockNames<-(data$cgenes$Block);
	
	#using the first array and block as the reference, and then normalize the rest array and block
	#array effects
	#+fcts<-matrix(rep(fct.list[["array"]][arrayNames],dim(data$C)[1]),
	#+			nrow=dim(data$C)[1], ncol=dim(data$C)[2],byrow=T);
	#+ndata$C<-ndata$C-fcts;
	
	#block effects
	#+fcts<-matrix(rep(fct.list[["block"]][blockNames],dim(data$C)[2]),
	#+			nrow=dim(data$C)[1], ncol=dim(data$C)[2],byrow=F);
	#+ndata$C<-ndata$C-fcts;
	#Done	
	ndata
}
#This is the function to do the normalization based on the RLM factor effects
#data is the data elist, assuming it has been log transformed 
#factor.effects is the data list, holding the parameters estimated by RLM, 
# this is the output of readFactors function.
# 
#coding: Simple or Dummy or Deviation
#Dummy coding, means using the first group as the reference, all 
#things are normalized towards this group.
#'@export
sub.effects.Dummy<-function(ndata,factor.effects)
{
	#check the input 
	if(missing(ndata))
	{
		stop("***Error: input data missing, please specify!!");
	}
	if(class(ndata)!="EListRaw")
	{
		stop("***Error: unknown input data format, please use EListRaw data"); 
	}
	#check the input 
	if(missing(factor.effects))
	{
		stop("***Error: input data missing, please specify!!");
	}
	if(class(factor.effects)!="list")
	{
		stop("***Error: unknown input data format, please use data list"); 
	}

	#----first doing normalization of E field (features)
	arrayNames<-colnames(ndata$E);
	blockNames<-(ndata$genes$Block);
	
	#using the first array and block as the reference, and then normalize the rest array and block
	#array effects
	fcts<-matrix(rep(factor.effects[["array"]][arrayNames],dim(ndata$E)[1]),
				nrow=dim(ndata$E)[1], ncol=dim(ndata$E)[2],byrow=T);
	ndata$E<-ndata$E-fcts;
	
	#block effects
	fcts<-matrix(rep(factor.effects[["block"]][blockNames],dim(ndata$E)[2]),
				nrow=dim(ndata$E)[1], ncol=dim(ndata$E)[2],byrow=F);
	ndata$E<-ndata$E-fcts;
	
	#-----now doing the control too
	arrayNames<-colnames(ndata$C);
	blockNames<-(ndata$cgenes$Block);
	
	#using the first array and block as the reference, and then normalize the rest array and block
	#array effects
	fcts<-matrix(rep(factor.effects[["array"]][arrayNames],dim(ndata$C)[1]),
				nrow=dim(ndata$C)[1], ncol=dim(ndata$C)[2],byrow=T);
	ndata$C<-ndata$C-fcts;
	
	#block effects
	fcts<-matrix(rep(factor.effects[["block"]][blockNames],dim(ndata$C)[2]),
				nrow=dim(ndata$C)[1], ncol=dim(ndata$C)[2],byrow=F);
	ndata$C<-ndata$C-fcts;
	#Done	
	ndata

} 
#
#data is the data elist, assuming it has been --->log transformed<---- 
#factor.effects is the data list, holding the parameters estimated by RLM, 
# this is the output of readFactors function. 
#coding: Simple or Dummy or Deviation
#Simple coding: means we are using the Grand mean as the intercept, everything is towards the
#			the reference (the first group). 
#'@export
sub.effects.Simple<-function(ndata,factor.effects)
{
	#check the input 
	if(missing(ndata))
	{
		stop("***Error: input data missing, please specify!!");
	}
	if(class(ndata)!="EListRaw")
	{
		stop("***Error: unknown input data format, please use EListRaw data"); 
	}
	#check the input 
	if(missing(factor.effects))
	{
		stop("***Error: input data missing, please specify!!");
	}
	if(class(factor.effects)!="list")
	{
		stop("***Error: unknown input data format, please use data list"); 
	}
	#cat("a1\n");
	#----first doing normalization of E field (features)
	arrayNames<-colnames(ndata$E);
	blockNames<-as.character(ndata$genes$Block);
	
	#using grand mean as the reference, and then normalize the rest array and block
	#array effects
	fct<-factor.effects[["array"]][-1] #get rid of zero one in the begin, the rest ones are involved in the calculation
	#need to make the coding matrix
	array.simple<-matrix(1/length(unique(arrayNames)), nrow=length(unique(arrayNames)),ncol=length(unique(arrayNames))-1,byrow=TRUE);
	###	
	array.simple<-contr.treatment(length(unique(arrayNames)))-array.simple;
	array.simple<-array.simple*matrix(rep(fct,dim(array.simple)[1]),nrow=dim(array.simple)[1],ncol=dim(array.simple)[2],byrow=TRUE)
	eff.array<-apply(array.simple,1, sum)  #->>need to add the grand mean, without addition, it still OK, but shifted 
	names(eff.array)<-unique(arrayNames);
	fcts<-matrix(rep(eff.array[arrayNames],dim(ndata$E)[1]),
				nrow=dim(ndata$E)[1], ncol=dim(ndata$E)[2],byrow=T);
	ndata$E<-ndata$E-fcts;
	#cat("a2\n");
	#------block effects
	block.simple<-matrix(1/length(unique(blockNames)), nrow=length(unique(blockNames)),ncol=length(unique(blockNames))-1,byrow=TRUE);
	fct<-factor.effects[["block"]][-1] #get rid of zero one in the begin, the rest ones are involved in the calculation

	block.simple<-contr.treatment(length(unique(blockNames)))-block.simple;
	block.simple<-block.simple*matrix(rep(fct,dim(block.simple)[1]),nrow=dim(block.simple)[1],ncol=dim(block.simple)-1,byrow=TRUE)
	eff.block<-apply(block.simple,1, sum)  #->>>add grand mean, without addition, it still OK, but shifted 
	names(eff.block)<-unique(blockNames)
	fcts<-matrix(rep(eff.block[blockNames],dim(ndata$E)[2]),
				nrow=dim(ndata$E)[1], ncol=dim(ndata$E)[2],byrow=F);
	ndata$E<-ndata$E-fcts;
	#cat("a3\n");
	#***********-----now doing the control too
	arrayNames<-colnames(ndata$C);
	blockNames<-as.character(ndata$cgenes$Block);;
	
	#array effects
	fct<-factor.effects[["array"]][-1] #get rid of zero one in the begin, the rest ones are involved in the calculation
	#need to make the coding matrix
	array.simple<-matrix(1/length(unique(arrayNames)), nrow=length(unique(arrayNames)),ncol=length(unique(arrayNames))-1,byrow=TRUE);
	###	
	array.simple<-contr.treatment(length(unique(arrayNames)))-array.simple;
	array.simple<-array.simple*matrix(rep(fct,dim(array.simple)[1]),nrow=dim(array.simple)[1],ncol=dim(array.simple)[2],byrow=TRUE)
	eff.array<-apply(array.simple,1, sum)
	names(eff.array)<-unique(arrayNames);
	fcts<-matrix(rep(eff.array[arrayNames],dim(ndata$C)[1]),
				nrow=dim(ndata$C)[1], ncol=dim(ndata$C)[2],byrow=T);
	ndata$C<-ndata$C-fcts;
	#cat("a4\n");
	#------block effects
	block.simple<-matrix(1/length(unique(blockNames)), nrow=length(unique(blockNames)),ncol=length(unique(blockNames))-1,byrow=TRUE);
	fct<-factor.effects[["block"]][-1] #get rid of zero one in the begin, the rest ones are involved in the calculation

	block.simple<-contr.treatment(length(unique(blockNames)))-block.simple;
	block.simple<-block.simple*matrix(rep(fct,dim(block.simple)[1]),nrow=dim(block.simple)[1],ncol=dim(block.simple)-1,byrow=TRUE)
	eff.block<-apply(block.simple,1, sum)
	names(eff.block)<-unique(blockNames)
	fcts<-matrix(rep(eff.block[blockNames],dim(ndata$C)[2]),
				nrow=dim(ndata$C)[1], ncol=dim(ndata$C)[2],byrow=F);
	ndata$C<- ndata$C-fcts;
	#cat("a5\n");
	#Done	
	ndata
}

#
#data is the data elist, assuming it has been --->log transformed<---- 
#factor.effects is the data list, holding the parameters estimated by RLM,
# array, block, intercept 
# this is the output of readFactors function. 
#coding: Simple or Dummy or Deviation？mainly Deviation(???, at least for sure doing deviation.)
#Deviation coding: means we are using the Grand cell mean as the reference, everything is towards the
#			the grand mean. 
#Simple coding: similar to Dummy coding, still using the first factor group as the reference, and all other 
#		groups compare to this reference. the only different is the regression intercept, which is the grand cell mean
#		instead of the reference cell mean.
#'@export
sub.effects<-function(ndata,factor.effects, coding=c("Deviation", "Simple","Dummy"))
{
	#check the input 
	if(missing(ndata))
	{
		stop("***Error: input data missing, please specify!!");
	}
	if(class(ndata)!="EListRaw")
	{
		stop("***Error: unknown input data format, please use EListRaw data"); 
	}
	#check the input 
	if(missing(factor.effects))
	{
		stop("***Error: input data missing, please specify!!");
	}
	if(class(factor.effects)!="list")
	{
		stop("***Error: unknown input data format, please use data list"); 
	}
	coding<-match.arg(coding);
	
	#cat("a1\n");
	#----first doing normalization of E field (features)
	arrayNames<-colnames(ndata$E);
	blockNames<-as.character(ndata$genes$Block);
	
	#using grand mean as the reference, and then normalize the rest array and block
	#array effects
	fct<-factor.effects[["array"]] #get rid of zero one in the begin, the rest ones are involved in the calculation
	#need to make the coding matrix
	array.mtx<-factor.effects[["contrasts"]][["array"]];
	rnames<-rownames(array.mtx);
	###	
	#array.simple<-contr.treatment(length(unique(arrayNames)))-array.simple;
	array.mtx<-array.mtx*matrix(rep(fct,dim(array.mtx)[1]),nrow=dim(array.mtx)[1],ncol=dim(array.mtx)[2],byrow=TRUE);
	eff.array<-apply(array.mtx,1, sum)  #->>need to add the grand mean, without addition, it still OK, but shifted 
	#if(coding=="Simple")
	#{
	#	eff.array<-eff.array+rep(factor.effects[["intercept"]],length(eff.array));
	#}
	names(eff.array)<-rnames
	fcts<-matrix(rep(eff.array[arrayNames],dim(ndata$E)[1]),
				nrow=dim(ndata$E)[1], ncol=dim(ndata$E)[2],byrow=T);
	ndata$E<-ndata$E-fcts;
	#cat("a2\n");

	#------block effects
	block.mtx<-factor.effects[["contrasts"]][["block"]];
    rnames<-rownames(block.mtx);
	fct<-factor.effects[["block"]] #get rid of zero one in the begin, the rest ones are involved in the calculation

	#block.mtx<-contr.treatment(length(unique(blockNames)))-block.simple;
	block.mtx<-block.mtx*matrix(rep(fct,dim(block.mtx)[1]),nrow=dim(block.mtx)[1],ncol=dim(block.mtx)-1,byrow=TRUE)
	eff.block<-apply(block.mtx,1, sum)  #->>>add grand mean, without addition, it still OK, but shifted 
	
	#if(coding=="Simple")
	#{
	#	eff.block<-eff.block+rep(factor.effects[["intercept"]],length(eff.array));
	#}
	names(eff.block)<-rnames;#unique(blockNames)
	fcts<-matrix(rep(eff.block[blockNames],dim(ndata$E)[2]),
				nrow=dim(ndata$E)[1], ncol=dim(ndata$E)[2],byrow=F);
	ndata$E<-ndata$E-fcts;
	if(coding=="Simple")
	{
		ndata$E<-ndata$E+matrix(factor.effects[["intercept"]], nrow=dim(ndata$E)[1],ncol=dim(ndata$E)[2]);
			#add back the grand mean (in Simple coding the intercept is the grand mean, otheriwse it is shift towards too low.
	}
	#cat("a3\n");
	
	
	#***********-----now doing the control too
	arrayNames<-colnames(ndata$C);
	blockNames<-as.character(ndata$cgenes$Block);
	
	#array effects
	#note the eff.array were calculated in the above E fields
	fcts<-matrix(rep(eff.array[arrayNames],dim(ndata$C)[1]),
				nrow=dim(ndata$C)[1], ncol=dim(ndata$C)[2],byrow=T);
	ndata$C<-ndata$C-fcts;
	#cat("a4\n");
	#------block effects
	fcts<-matrix(rep(eff.block[blockNames],dim(ndata$C)[2]),
				nrow=dim(ndata$C)[1], ncol=dim(ndata$C)[2],byrow=F);
	ndata$C<- ndata$C-fcts;
	
	if(coding=="Simple")
	{
		ndata$C<-ndata$C+matrix(factor.effects[["intercept"]], nrow=dim(ndata$C)[1],ncol=dim(ndata$C)[2]);
			#add back the grand mean (in Simple coding the intercept is the grand mean, otheriwse it is shift towards too low.
	}
	#cat("a5\n");
	#Done	
	ndata
} 
 

#parseParameters
#'@title Reading the coefficients from RLM fit
#'
#'@description read the estimated coefficients from RLM fits and parse them into the factor effects.
#'
#'@details Read the RLM fit object and parse the coefficients from the fit according to the effects. 
#' Arrange them into factor effects. Currently only (first) two factors are parsed, 
#' assuming they are "array" and "block". Leave "feature" factors unread even if they are 
#' estimated. When parsing, we don't anything but simply parse the values into different factors.
#' Of course, We do name the 
#' output factor levels. (Well, it seems likely it doesn't matter how 
#'	to call them. Only the order matters!!!need to confirm.)
#'
#'@param fit rlm object the input RLM fit object of class rlm and lm
#'@param factors string vector holding the names of the factors/variables used to fit
#'	RLM. It could be either c("array", "block") or c("array", "block", "feature"). 
#'	The former case is when there is only one control protein to fit RLM, and the 
#'	latter there are two or more controls selected for fitting.
#'@param factor.levels the numeric vector used to indicate how many levels there are
#'	in each factor. It should have the same length and order as those of factors vector.
# # ->>not using now, obseleted '@param factor.first the string vector that hold the name string of the first 
# #'	level of each factor. The reason we need these name strings is that when 
# #'	fitting RLM with a categorical variable, the coefficients are only estimated
## '	for other levels relative to the first one. The information of first one is 
## '	contained in the intercept by default. So the names of levels of the factors
# #'	excluding the first levels can be extracted from coefficients, but we need
# #'	the names of the first levels to be provided. While in the case of "Deviation"
# #'	and "Simple" coding, the intercept contains the grand mean instead of the first level. 
#'
#'@return a list containing the named numeric vectors holding the 
#'	factor effects estimated by RLM. Currently, only two factors,
#'	"array" and "block", are extracted and returned. Each vector
#' contained all factor estimator (the number of factor estiamtor is one less
#'	than the factor levels). We will also return the intercept, since for "Simple"
#'	and "Deviation" coding, this is the grand mean and we need it to calculate
#'	the factor effects!!!.
#'
# '@example
#'
#'@export
 
readFactors<-function(fit, factors, factor.levels #, factor.first
				)
{
	if(missing(fit))
	{
		stop("***Error: input fitting object missing, please specify!!");
	}
	if(all((class(fit)!=c("rlm","lm"))))
	{
		stop("***Error: unknown input object format, please use rlm fitting object"); 
	}
	if(missing(factor.levels))
	{
		stop("***Error: input factor.levels missing, please specify!!");
	}
	if(class(factor.levels)!="integer"&&class(factor.levels)!="numeric")
	{
		stop("***Error: unknown input object format, please use integer vector"); 
	}
	if(missing(factors))
	{
		stop("***Error: input factors missing, please specify!!");
	}
	if(class(factors)!="character")
	{
		stop("***Error: unknown input factors format, please use character vectors"); 
	}
	#if(missing(factor.first))
	#{
	#	stop("***Error: input factor.first missing, please specify!!");
	#}
	#if(class(factor.first)!="character")
	#{
	#	stop("***Error: unknown input factor.first format, please use character vectors"); 
	#}	
	#coding<-match.args
	#now start parsing
	if(length(factors)!=length(factor.levels))#||length(factor.levels)!=length(factor.first))
	{
		stop("**ERROR: the input factors and factor.levels and factor.first vectors are not of equal length, please check!!!");
	}
	
	cff<-fit$coefficients
	if(length(cff) < (sum(factor.levels)-length(factor.levels)))
	{
		stop("ERROR: number of coefficients by RLM fitting does not aggress with the input factor levels, please check!!"); 
	}
	cff<-cff[-1]; #the first one is the intercept, the base level.
	ret.list<-vector("list",length=length(factors)+2)
	#coding.matrix<-vector("list",length=length(factors))
	for(i in 1:length(factors))###assuming the order of factors is the same as the one for fitted parameters
	{
			#getting array first
			fct_prefix<-factors[i];
			fct.values<-cff[1:(factor.levels[i]-1)];
			fct.names<-names(cff)[1:(factor.levels[i]-1)];
			fct.names<-substr(fct.names,start=nchar(fct_prefix)+1,stop=100000L);
			#fct.values<-c(0,fct.values)
			#fct.names<-c(factor.first[i],fct.names)
			names(fct.values)<-fct.names;
			
			cff<-cff[-(1:(factor.levels[i]-1))]; #remove the ones that have been parsed, get ready for the next round
			ret.list[[i]]<-fct.values;
			#coding.matrix[[i]]<-contrasts(
	}
	ret.list[[length(factors)+1]]<-fit$coefficients[1];
	ret.list[[length(factors)+2]]<-fit$contrasts;
	names(ret.list)<-c(factors,"intercept","contrasts");
	
	ret.list
}
###this is the function to format the elist object into
###data.frame that can be easily linearly fitted.
###the model is like:
###  Y_ijkm=a_i+b_j+f_k+epsilon_ijkm
###			a: array; b: block; f: feature; epsilon: Gaussian errors
##in this function, we need to convert the control protein data 
## matrix into data frame, so to do linear fitting
#'@title control protein data to data frame for RLM fitting
#'@description the function used to convert the selected control protein
#'	matrix into a data frame that is ready for RLM fitting
#'@details It takes the control protein data matrix and control gene description
#'	table and then select a subset of control protein in order to fit RLM to estimate
#'	array and block effects according to the linear model based on work by Sboner et al.
#'	Briefly, the model is like \cr 
#'	\verb{   } Y_ijkm=a_i+b_j+f_k+epsilon_ijkm
#'			\cr a: array; b: block; f: feature; epsilon: Gaussian errors
#'@param cmatrix data matrix containing the target/feature protein expression levels
#'	arranged in gene (row) by array (column) format. It is the "C" field of the protoArray
#'	EListRaw object. 
#'@param cgene data frame containing the meta data about the control proteins. 
#'	It is the "cgenes" field of field of the protoArray	EListRaw object. 
#'@param control string the names to select the control genes for fitting. It could 
#'	be one single string or a regexp term to select a group of control proteins.
#'
#'@param log logic used to indicate whether logarithm transformation of data is
#'	necessary  
#'@param log.base numeric the base for the logarithmic transformation if necessary.
#'@return a data frame containing the data and variables for RLM fitting. The 
#'	columns named in a fixed structure. \cr
#'	\tabular{ccccc}{
#'		\tab ----------- \tab ----------- \tab ----------- \tab ----------- \cr 
#' 		\tab Y \tab array \tab block \tab feature \cr 
#'		\tab ----------- \tab ----------- \tab ----------- \tab ----------- \cr  
#'	}
#'	The last column is optional depending on the data input.
## '@example  
#'@export

cmatrix2dataframe<-function(cmatrix, cgene, control, log=TRUE, log.base=exp(1))
{
	#dtf.C<-data.frame()
	#check the input 
	if(missing(cmatrix))
	{
		stop("***Error: input data missing, please specify!!");
	}
	if(class(cmatrix)!="matrix")
	{
		stop("***Error: unknown input data format, please use data matrix"); 
	}
	if(missing(cgene))
	{
		stop("***Error: input cgene table missing, please specify!!");
	}
	if(class(cgene)!="data.frame")
	{
		stop("***Error: unknown input cgene format, please use data frame"); 
	}
	#now get a, array effects and feature data frame first, block effects
	#first get the data set ready
	if(!missing(control)&&!is.null(control))
	{
		idx<-grep(control,cgene$Name);
		if(length(idx)==0)
		{
			stop(paste("Error: the specified gene name(s) \"", control, "\" returns no entries, please check input!!!",sep=""));
		}
		cgene<-cgene[idx,];
		cmatrix<-cmatrix[idx,];
	}
	#first get unique feature set
	f.names<-unique(cgene$Name);
	
	#get array names
	array.names<-colnames(cmatrix)
	if(is.null(array.names))
	{
		array.names<-paste("a",1:dim(cmatrix)[2],sep="_");
	}
	#array.name<-paste(array.name,i,sep="_");
	#array.name_i<-rep(array.name,dim(data$E)[1]);
		
	for(i in 1:length(f.names))
	{	
		#get each gene for all
		cName<-f.names[i];
		#get all this ctrl protein out
		cfeature.indices<-which(cgene$Name==cName)
		
		#put the things into the data frame
		curr_dtf<-NULL;
		for(j in 1:length(cfeature.indices))
		{
			curr_dtf<-data.frame("Y"=cmatrix[cfeature.indices[j],],"feature"=cName,
							"array"=array.names,"block"=cgene[cfeature.indices[j],"Block"])
			if(i==1&&j==1){
				dtf.C<-curr_dtf;
			} else {
				dtf.C<-rbind(dtf.C,curr_dtf);
			}
		}#end of for loop for each feature

	}##end of for loop for all features
	dtf.C$block<-as.factor(dtf.C$block)
	if(log)
	{
		dtf.C$Y<-log(dtf.C$Y, base=log.base)
	}
	dtf.C
}

#### function to run RLM fitting with a sboner model
###  Y_ijkm=a_i+b_j+f_k+epsilon_ijkm
###			a: array; b: block; f: feature; epsilon: Gaussian errors
#####in this function, we fit control features to estimate 
###, a_i and b_j, then 
## model, no interaction
## also all 3 variables are categorical
## need also to check when there is one feature, ie only one control protein selected!!
##			 
#'@export
#
fit.RLM<-function(data, coding=c("Deviation","Dummy","Simple"))
{
	if(missing(data))
	{
		stop("***Error: input data missing, please specify!!");
	}
	if(class(data)!="data.frame")
	{
		stop("***Error: unknown input data format, please use data frame"); 
	}
	coding<-match.arg(coding);
	
	
	#check for data consistency
	fieldNames<-names(data);
	if(!all(is.element(c("array", "block","feature"),fieldNames)))
	{
		stop("Error: some fields are not in the input data frame,please check the input data"); 
	}
	#contrasts(data$array)<-contr.treatment(length(unique(data$array)));
	#contrasts(data$block)<-contr.treatment(length(unique(data$block)));
	#if we want to do simply, we need to construct the matrix
	#if we want to do dummy, just don't do anything, since it is the default
	#####***** if we want to call explicitly to do dummy coding, we have to set the coding matrix colnames, 
	#####*****otherwise it is replaced with variable name suffixed with number , eg array2. 
		
	if(coding=="Simple"){ #the default is dummy
		array.simple<-matrix(1/length(unique(data$array)), nrow=length(unique(data$array)),ncol=length(unique(data$array))-1,byrow=TRUE);
		block.simple<-matrix(1/length(unique(data$block)), nrow=length(unique(data$block)),ncol=length(unique(data$block))-1,byrow=TRUE);
		feature.simple<-NULL;
		if(length(unique(data$feature))!=1)
		{
			feature.simple<-matrix(1/length(unique(data$feature)), nrow=length(unique(data$feature)),ncol=length(unique(data$feature))-1,byrow=TRUE);
		}
		
		c<-contrasts(data$array)-array.simple;
		colnames(c)<-unique(data$array)[-1];
		contrasts(data$array)<-c;
		
		c<-contrasts(data$block)-block.simple;
		colnames(c)<-unique(data$block)[-1];
		contrasts(data$block)<-c;
		
		if(length(unique(data$feature))!=1)
		{
			c<-contrasts(data$feature)-feature.simple;
			colnames(c)<-unique(data$feature)[-1];
			contrasts(data$feature)<-c;
		}
	}
	
	#this following is for "Deviation" contrasts
	if(coding=="Deviation"){ #the default is dummy
		array.deviation<-contr.sum(length(unique(data$array)));
		block.deviation<-contr.sum(length(unique(data$block)));# nrow=length(unique(data$block)),ncol=length(unique(data$block))-1,byrow=TRUE);
				
		if(length(unique(data$feature))!=1)
		{
			feature.deviation<-contr.sum(length(unique(data$feature)));
		}
		c<-array.deviation;
		colnames(c)<-unique(data$array)[-1];
		contrasts(data$array)<-c;
		
		c<-block.deviation;
		colnames(c)<-unique(data$block)[-1];
		contrasts(data$block)<-c;
				
		if(length(unique(data$feature))!=1)
		{
			c<-feature.deviation;
			colnames(c)<-unique(data$feature)[-1];
			contrasts(data$feature)<-c

		}
	}
	#contrasts(data$feature)<-contr.treatment(length(unique(data$feature)));
	fit.rlm<-NULL;
	cat(paste("Start doing the RLM fitting -- ",coding," coding:\n",sep=""));
	if(length(unique(data$feature))==1)
	{
		cat("\tRLM model: only one control protein selected for fitting!\n") ;
		fit.rlm<-rlm(Y~array+block,data=data,psi=psi.huber ,maxit=100);
	} else {
		fit.rlm<-rlm(Y~array+block+feature,data=data,psi=psi.huber,maxit=100);
	}
	cat("Done!!!\n");
	return(fit.rlm)
}

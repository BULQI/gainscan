##########################################
##the code here is used to run simulation tests for 
##checking the model and methods. 
##1.testing whether the simulated data can be used to
##	detecting the nonZero interaction effects. ----it should
##2.checking for the heterscedasticity cases vs. homoscedasticity
##3.checking for validity of the modeling without control data
##4.all above first using the smaller number of genes---- to keep 
##	the memory usage low so as to using the R built-in linear 
##	function for testing
##
##	Oct 1st, 2016 Feng@BU
##	see ModelTesting_run1.R for equal variance cases
##
##
##@Oct 23, 2016 by Feng
##ModelTesting code file #3
##simulate datafor unequal variance cases
##here, we need to do the transformation!!!
##since clearly the unequal variances have lower true discovery rate.
## To transform we need to calculate the expected variance using bayesian approach
## (smyth 2004), then transform to stabalize the variance by 
##	 var_hat=(d0S0^2+diSi^2)/(d0+di)
##	 
##		Y_t=(Y-Y_bar)/sqrt(var_hat)+Y_bar
##by doing this, we kind of stablized the variance, but did not 
##changed the variance
##
##
##in the end we do the linear regression 
##check the assumption
##as well as the evaluate the performance
##
##########################################
library(ARPPA)
##first checking for homoscedastical dataset
#setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_NegativeCon_Trans")
setwd("E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_NegativeCon_Trans")
nGene<-2000
nTreatment<-2   #number of different beta
sampleSize<-5 ##this is the number of repeats for each group
alpha.mean<-0  #variance for alpha prior
alpha.sigma<-3
beta.mean<-0
beta.sigma<-2   #variance for beta prior

gamma.sigma<-10  # the unscaled factor for gamma given gamma<=>0
prob.nonZero<-0.02
#priors for variance distribution
d0<-5
s0<-2
repeats<-100

#this time of around, we need a fixed variance for the gaussian variance
#epsilon.si<-2 #not use for this case
#call it
set.seed(2004);
#simulate *****the equal variance and !!!!!!!negative control data
#E-equal variance; N-negative control
i<-1
for(i in c(1:repeats))
{
cat("doing ", i, "/",repeats," analyses....\n")
	dataExp_EN_List<-simulateExpression(nGene, nTreatment, sampleSize,
				control.negative=TRUE, control.isotype=FALSE,
				control.index=c(1),
				alpha.mean=alpha.mean, alpha.sigma=alpha.sigma,
				beta.mean=beta.mean, beta.sigma=beta.sigma,
				prob.nonZero=prob.nonZero, gamma.sigma=gamma.sigma,
				#epsilon.si=epsilon.si,
				epsilon.d0=d0, epsilon.s0=s0
				)
#calculate the sample variance from the data
dataExp_EN<-dataExp_EN_List[[1]]
svList<-sampleVariance(dataExp_EN,nTreatment, sampleSize)
#reformat it into a vector
s2<-as.vector(as.matrix(svList))
df<-as.vector(as.matrix(sampleSizes(dataExp_EN,nTreatment, sampleSize)))-1
r2<-sampleVariancePrior(s2,df)

#now do the transformation to stabilize the variance
dataTransformed<-transformData(dataExp_EN,nTreatment, sampleSize)
#we have the data, what to do.
#reformat the data from data matrix to 
#dataframe for the linear regress
dExp_EN<-matrix2dframe(dataTransformed, nTreatment, sampleSize);#transformed
#dExp_EN<-matrix2dframe(dataExp_EN, nTreatment, sampleSize); #untransformed 


#now do the linear regression with interaction 

	lreg_EN<-lm(exp~gene*group, data=dExp_EN) 
	lregSm<-summary(lreg_EN)
	lstToSave<-list("data"=dataExp_EN_List, "lregCoef"=lregSm$coefficients, "var.prior"=r2);
	save(lstToSave, file=paste("result_",i,".RData",sep=""));

}

###############
set.seed(2005)
setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_IsotypeCon_Trans")
j<-1
for(j in c(1:repeats))
{
cat("doing ", j, "/",repeats," analyses....\n")
#simulate *****the equal variance and !!!!!!!isotype control data
#E-equal; I-isotype

dataExp_EI_List<-simulateExpression(nGene, nTreatment, sampleSize,
				control.negative=FALSE, control.isotype=TRUE,
				control.index=c(1),
				alpha.mean=alpha.mean, alpha.sigma=alpha.sigma,
				beta.mean=beta.mean, beta.sigma=beta.sigma,
				prob.nonZero=prob.nonZero, gamma.sigma=gamma.sigma,
				#epsilon.si=epsilon.si #
				epsilon.d0=d0, epsilon.s0=s0
				)
dataExp_EI<-dataExp_EI_List[[1]]
#calculate the sample variance from the data
svList<-sampleVariance(dataExp_EI,nTreatment, sampleSize)
#reformat it into a vector
s2<-as.vector(as.matrix(svList))
df<-as.vector(as.matrix(sampleSizes(dataExp_EI,nTreatment, sampleSize)))-1
r2<-sampleVariancePrior(s2,df)

#now do the transformation
dataTransformed<-transformData(dataExp_EI,nTreatment, sampleSize)

#we have the data, what to do.
#reformat the data from data matrix to 
#dataframe for the linear regress
dExp_EI<-matrix2dframe(dataTransformed, nTreatment, sampleSize);

#now do the linear regression with interaction
if(j>65)
 {  
	lreg_EI<-lm(exp~gene*group, data=dExp_EI)
	lregSm_EI<-summary(lreg_EI)
	lstToSave_EI<-list("data"=dataExp_EI_List, "lregCoef"=lregSm_EI$coefficients, "var.prior"=r2);
	save(lstToSave_EI, file=paste("result_",j,".RData",sep=""));
	cat("Done!!\n")
 }
}

#####################
set.seed(2006)
setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_UnControl_Trans")
k<-1
for(k in c(1:repeats))
{
cat("doing ", k, "/",repeats," analyses....\n")
#simulate *****the equal variance and !!!!!!!isotype control data
#E-equal; NC-NC

dataExp_ENC_List<-simulateExpression(nGene, nTreatment, sampleSize,
				control.negative=FALSE, control.isotype=FALSE,
				#control.index=c(1),
				alpha.mean=alpha.mean, alpha.sigma=alpha.sigma,
				beta.mean=beta.mean, beta.sigma=beta.sigma,
				prob.nonZero=prob.nonZero, gamma.sigma=gamma.sigma,
				#epsilon.si=epsilon.si 
				epsilon.d0=d0, epsilon.s0=s0
				)
dataExp_ENC<-dataExp_ENC_List[[1]]
#calculate the sample variance from the data
svList<-sampleVariance(dataExp_ENC,nTreatment, sampleSize)
#reformat it into a vector
s2<-as.vector(as.matrix(svList))
df<-as.vector(as.matrix(sampleSizes(dataExp_ENC,nTreatment, sampleSize)))-1
r2<-sampleVariancePrior(s2,df)

#now do the transformation
dataTransformed<-transformData(dataExp_ENC,nTreatment, sampleSize)

#we have the data, what to do.
#reformat the data from data matrix to 
#dataframe for the linear regress
dExp_ENC<-matrix2dframe(dataTransformed, nTreatment, sampleSize);

#now do the linear regression with interaction 
if(k>64)
	{
	lreg_ENC<-lm(exp~gene*group, data=dExp_ENC)
	lregSm_ENC<-summary(lreg_ENC)
	lstToSave_ENC<-list("data"=dataExp_ENC_List, "lregCoef"=lregSm_ENC$coefficients, "var.prior"=r2);
	save(lstToSave_ENC, file=paste("result_",k,".RData",sep=""));
	cat("Done!!\n")
	}
}
#now summary the data first
aggregate(exp~gene*group, dExp_ENC, FUN=mean)

###dianostic plots
op<-par(mfrow=c(2,2))
plot(lreg_EN)
par(op)

##############################################
#-------------------------------------------
##############################################
###the code below has not been used######
#===============================================================

#================data analysis=======
###start doing the summary statistics of the repeated runs
#
#loading the data
repeats<-100
setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_NegativeCon");
#setwd("~/Desktop/arppr/EqualVar_NegativeCon/")
i<-1
portions<-seq(0.001, 0.5, 0.001)
stats_df<-data.frame("portions"=portions);
for(i in c(1:repeats))
{
	cat("reading the ",i, "/", repeats," data sets.....\n")
	load(paste("result_", i,".RData",sep=""))
	cat("count the correct ones.....\n")
	#depending on which data,
	lstData<-lstToSave
	rm(lstToSave)
	
	#try to know which one is nonZero as the true values
	gammaTrue<-lstData$data$gamma
	gammaTrueIndex<-which(gammaTrue[,2]!=0)
	
	#now need to sort the p-value of the 
	#and pick top n genes to be the significant ones as the analysis results
	#
	lr_coe<-lstData$lregCoef
	lr_coe_names<-rownames(lr_coe)
	lr_coe_gamma_index<-which(grepl("gene\\d*:group2",lr_coe_names));
	lr_coe_int<-lr_coe[lr_coe_gamma_index,]
	lr_coe_int_names<-lr_coe_names[lr_coe_gamma_index];

	#sort the array according to the prob
	lr_coe_int_sort_order<-order(lr_coe_int[,4]) #NOTE:the names is one more than its index, eg. 811 is gene812:group2
	#lr_coe_int_sort<-lr_coe_int[lr_coe_int_sort_order,]
	

	#now we have everything, just need to collect statistics
	p_correct<-portions;#initialize the vector
	for(j in c(1:length(portions)))
	{
		#for each portion, we need to check what is percentage to be correct
		numToPick<-floor(lstData$data$params[1]*portions[j])
		#picking from the order array
		genesToPick<-lr_coe_int_sort_order[c(1:numToPick)]+1

		#now we got the genes, just need to check whether the genes so far are the TRUE ones
		numOfCorrect<-sum(is.element(genesToPick, gammaTrueIndex));
		p_correct[j]<-numOfCorrect/numToPick
	}
	stats_df[,paste("repeat_",i)]<-p_correct;
}
#statistics with the stats_df array
mean_correct<-portions
max_correct<-portions
min_correct<-portions
std_correct<-portions
m_stats_df<-as.matrix(stats_df)
for(k in c(1:length(portions)))
{
	mean_correct[k]<-mean(m_stats_df[k,c(-1)])
	max_correct[k]<-max(m_stats_df[k,c(-1)])
	min_correct[k]<-min(m_stats_df[k,c(-1)])
	std_correct[k]<-sqrt(var(m_stats_df[k,c(-1)]))
}

stats_df[,"mean"]<-mean_correct
stats_df[,"min"]<-min_correct
stats_df[,"max"]<-max_correct
stats_df[,"std"]<-std_correct

plot(c(0.001,0.51),c(0.02, 1.1), type="n", main="true discovery rate for protein selection", 
	xlab="portion of protein selected", ylab="portion of true discovery", log="xy")
lines(stats_df[,1], stats_df[,"mean"], col=1, lty=1, lwd=2)
lines(stats_df[,1], stats_df[,"min"], col="grey", lty=2, lwd=1)
lines(stats_df[,1], stats_df[,"max"], col="grey", lty=2, lwd=1)
lines(c(0.01, 0.01), c(0.1,1.2), col=2, lty=3)

#show stats
stats_df[c(1:30),c("portions", "mean", "min", "max", "std")]

########=================Isotype control data analysis
#loading the data
repeats<-100
setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_IsotypeCon");

#setwd("~/Desktop/arppr/EqualVar_NegativeCon/")
i<-1
portions<-seq(0.001, 0.5, 0.001)
stats_df<-data.frame("portions"=portions);
for(i in c(1:repeats))
{
	cat("reading the ",i, "/", repeats," data sets.....\n")
	load(paste("result_", i,".RData",sep=""))
	cat("count the correct ones.....\n")
	#depending on which data,
	lstData<-lstToSave_EI
	rm(lstToSave_EI)
	
	#try to know which one is nonZero as the true values
	gammaTrue<-lstData$data$gamma
	gammaTrueIndex<-which(gammaTrue[,2]!=0)
	
	#now need to sort the p-value of the 
	#and pick top n genes to be the significant ones as the analysis results
	#
	lr_coe<-lstData$lregCoef
	lr_coe_names<-rownames(lr_coe)
	lr_coe_gamma_index<-which(grepl("gene\\d*:group2",lr_coe_names));
	lr_coe_int<-lr_coe[lr_coe_gamma_index,]
	lr_coe_int_names<-lr_coe_names[lr_coe_gamma_index];

	#sort the array according to the prob
	lr_coe_int_sort_order<-order(lr_coe_int[,4]) #NOTE:the names is one more than its index, eg. 811 is gene812:group2
	#lr_coe_int_sort<-lr_coe_int[lr_coe_int_sort_order,]
	

	#now we have everything, just need to collect statistics
	p_correct<-portions;#initialize the vector
	for(j in c(1:length(portions)))
	{
		#for each portion, we need to check what is percentage to be correct
		numToPick<-floor(lstData$data$params[1]*portions[j])
		#picking from the order array
		genesToPick<-lr_coe_int_sort_order[c(1:numToPick)]+1

		#now we got the genes, just need to check whether the genes so far are the TRUE ones
		numOfCorrect<-sum(is.element(genesToPick, gammaTrueIndex));
		p_correct[j]<-numOfCorrect/numToPick
	}
	stats_df[,paste("repeat_",i)]<-p_correct;
}
#statistics with the stats_df array
mean_correct<-portions
max_correct<-portions
min_correct<-portions
std_correct<-portions
m_stats_df<-as.matrix(stats_df)
for(k in c(1:length(portions)))
{
	mean_correct[k]<-mean(m_stats_df[k,c(-1)])
	max_correct[k]<-max(m_stats_df[k,c(-1)])
	min_correct[k]<-min(m_stats_df[k,c(-1)])
	std_correct[k]<-sqrt(var(m_stats_df[k,c(-1)]))
}

stats_df[,"mean"]<-mean_correct
stats_df[,"min"]<-min_correct
stats_df[,"max"]<-max_correct
stats_df[,"std"]<-std_correct

plot(c(0.001,0.51),c(0.02, 1.1), type="n", main="true discovery rate for protein selection", 
	xlab="portion of protein selected", ylab="portion of true discovery", log="xy")
lines(stats_df[,1], stats_df[,"mean"], col=1, lty=1, lwd=2)
lines(stats_df[,1], stats_df[,"min"], col="grey", lty=2, lwd=1)
lines(stats_df[,1], stats_df[,"max"], col="grey", lty=2, lwd=1)
lines(c(0.01, 0.01), c(0.02,1.2), col=2, lty=3)

#show stats
stats_df[c(1:30),c("portions", "mean", "min", "max", "std")]
 
 
 #######====================No control++++++++++

#loading the data
repeats<-100
setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_UnControl");
#setwd("~/Desktop/arppr/EqualVar_NegativeCon/")
i<-1
portions<-seq(0.001, 0.5, 0.001)
stats_df<-data.frame("portions"=portions);
for(i in c(1:repeats))
{
	cat("reading the ",i, "/", repeats," data sets.....\n")
	load(paste("result_", i,".RData",sep=""))
	cat("count the correct ones.....\n")
	#depending on which data,
	lstData<-lstToSave_ENC
	rm(lstToSave_ENC)
	
	#try to know which one is nonZero as the true values
	gammaTrue<-lstData$data$gamma
	gammaTrueIndex<-which(gammaTrue[,2]!=0)
	
	#now need to sort the p-value of the 
	#and pick top n genes to be the significant ones as the analysis results
	#
	lr_coe<-lstData$lregCoef
	lr_coe_names<-rownames(lr_coe)
	lr_coe_gamma_index<-which(grepl("gene\\d*:group2",lr_coe_names));
	lr_coe_int<-lr_coe[lr_coe_gamma_index,]
	lr_coe_int_names<-lr_coe_names[lr_coe_gamma_index];

	#sort the array according to the prob
	lr_coe_int_sort_order<-order(lr_coe_int[,4]) #NOTE:the names is one more than its index, eg. 811 is gene812:group2
	#lr_coe_int_sort<-lr_coe_int[lr_coe_int_sort_order,]
	

	#now we have everything, just need to collect statistics
	p_correct<-portions;#initialize the vector
	for(j in c(1:length(portions)))
	{
		#for each portion, we need to check what is percentage to be correct
		numToPick<-floor(lstData$data$params[1]*portions[j])
		#picking from the order array
		genesToPick<-lr_coe_int_sort_order[c(1:numToPick)]+1

		#now we got the genes, just need to check whether the genes so far are the TRUE ones
		numOfCorrect<-sum(is.element(genesToPick, gammaTrueIndex));
		p_correct[j]<-numOfCorrect/numToPick
	}
	stats_df[,paste("repeat_",i)]<-p_correct;
}
#statistics with the stats_df array
mean_correct<-portions
max_correct<-portions
min_correct<-portions
std_correct<-portions
m_stats_df<-as.matrix(stats_df)
for(k in c(1:length(portions)))
{
	mean_correct[k]<-mean(m_stats_df[k,c(-1)])
	max_correct[k]<-max(m_stats_df[k,c(-1)])
	min_correct[k]<-min(m_stats_df[k,c(-1)])
	std_correct[k]<-sqrt(var(m_stats_df[k,c(-1)]))
}

stats_df[,"mean"]<-mean_correct
stats_df[,"min"]<-min_correct
stats_df[,"max"]<-max_correct
stats_df[,"std"]<-std_correct

plot(c(0.001,0.51),c(0.02, 1.1), type="n", main="true discovery rate for protein selection", 
	xlab="portion of protein selected", ylab="portion of true discovery", log="")
lines(stats_df[,1], stats_df[,"mean"], col=1, lty=1, lwd=2)
lines(stats_df[,1], stats_df[,"min"], col="grey", lty=2, lwd=1)
lines(stats_df[,1], stats_df[,"max"], col="grey", lty=2, lwd=1)
lines(c(0.01, 0.01), c(0.02,1.2), col=2, lty=3)

#show stats
stats_df[c(1:30),c("portions", "mean", "min", "max", "std")]
 
 
 ####------------------------------------------------------------------------------------
 #####the following part is also for the ---no control---
 ###this is special for no control cases, because in this case
 ###the differential ones could possibly from both groups
 ###if we only count one group, this will lower the true positive rate
 ###here we need to check both groups for the true positive ones
 
 #loading the data
 
 repeats<-100
setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_UnControl");
#setwd("~/Desktop/arppr/EqualVar_NegativeCon/")
i<-1
portions<-seq(0.001, 0.5, 0.001)
stats_df<-data.frame("portions"=portions);
for(i in c(1:repeats))
{
	cat("reading the ",i, "/", repeats," data sets.....\n")
	load(paste("result_", i,".RData",sep=""))
	cat("count the correct ones.....\n")
	#depending on which data,
	lstData<-lstToSave_ENC
	rm(lstToSave_ENC)
	
	#try to know which one is nonZero as the true values
	gammaTrue<-lstData$data$gamma
	gammaTrueIndex2<-which(gammaTrue[,2]!=0)
	gammaTrueIndex1<-which(gammaTrue[,1]!=0)
	#now need to sort the p-value of the 
	#and pick top n genes to be the significant ones as the analysis results
	#
	lr_coe<-lstData$lregCoef
	lr_coe_names<-rownames(lr_coe)
	lr_coe_gamma_index<-which(grepl("gene\\d*:group2",lr_coe_names));
	lr_coe_int<-lr_coe[lr_coe_gamma_index,]
	lr_coe_int_names<-lr_coe_names[lr_coe_gamma_index];

	#sort the array according to the prob
	lr_coe_int_sort_order<-order(lr_coe_int[,4]) #NOTE:the names is one more than its index, eg. 811 is gene812:group2
	#lr_coe_int_sort<-lr_coe_int[lr_coe_int_sort_order,]
	

	#now we have everything, just need to collect statistics
	p_correct<-portions;#initialize the vector
	for(j in c(1:length(portions)))
	{
		#for each portion, we need to check what is percentage to be correct
		numToPick<-floor(lstData$data$params[1]*portions[j]*2)
		#picking from the order array
		genesToPick<-lr_coe_int_sort_order[c(1:numToPick)]+1

		#now we got the genes, just need to check whether the genes so far are the TRUE ones
		numOfCorrect1<-sum(is.element(genesToPick, gammaTrueIndex1));
		numOfCorrect2<-sum(is.element(genesToPick, gammaTrueIndex2));
		p_correct[j]<-(numOfCorrect1+numOfCorrect2)/numToPick
	}
	stats_df[,paste("repeat_",i)]<-p_correct;
}
#statistics with the stats_df array
mean_correct<-portions
max_correct<-portions
min_correct<-portions
std_correct<-portions
m_stats_df<-as.matrix(stats_df)
for(k in c(1:length(portions)))
{
	mean_correct[k]<-mean(m_stats_df[k,c(-1)])
	max_correct[k]<-max(m_stats_df[k,c(-1)])
	min_correct[k]<-min(m_stats_df[k,c(-1)])
	std_correct[k]<-sqrt(var(m_stats_df[k,c(-1)]))
}

stats_df[,"mean"]<-mean_correct
stats_df[,"min"]<-min_correct
stats_df[,"max"]<-max_correct
stats_df[,"std"]<-std_correct

plot(c(0.001,0.51),c(0.02, 1.1), type="n", main="true discovery rate for protein selection", 
	xlab="portion of protein selected", ylab="portion of true discovery", log="xy")
lines(stats_df[,1], stats_df[,"mean"], col=1, lty=1, lwd=2)
lines(stats_df[,1], stats_df[,"min"], col="grey", lty=2, lwd=1)
lines(stats_df[,1], stats_df[,"max"], col="grey", lty=2, lwd=1)
lines(c(0.01, 0.01), c(0.02,1.2), col=2, lty=3)

#show stats
stats_df[c(1:30),c("portions", "mean", "min", "max", "std")]

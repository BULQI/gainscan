##########################################
## Dec 2016 code
##the code here are used to run simulation tests for 
##checking the model and methods, especially check the assumptions. 
##Read/load previous saved data including the linear regression
## and then plot for diagnostics
##
## =========> here we only use the negative control cases as the example to show the assumptions
##   feng@BU
##
##   ------------------
##	Dec 16, 2916, update
##  for the data analysis for performance, we changed to use the function for analysis
##		see the functionAnalysis.r for function definition.
##  we also save the original code to ModelTesting_run4.R.old
##########################################
library(ARPPA)
##first checking for weighted unEqualVar dataset -- Negative control
setwd("E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_NegConWt")
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
	sMean<-sampleMean(dataExp_EN, nTreatment, sampleSize)
	sVar<-sampleVariance(dataExp_EN, nTreatment, sampleSize)
	sdf<-sampleSizes(dataExp_EN, nTreatment, sampleSize)-1
	#get the sample variance prior by assuming a scale-chisquare
	s2<-as.vector(as.matrix(sVar))
	df<-as.vector(as.matrix(sdf))
	
	r2<-sampleVariancePrior(s2,df)
	s02<-r2[[2]]
	df0<-r2[[1]]
	
	#calculate the posterior variance
	var_pos<-sampleVariancePosterior(sdf,sVar,  df0, s02 )
	var_pos_vec<-(as.vector(t(var_pos)))
	var_pos_vec_all<-rep(var_pos_vec,rep(5, length(var_pos_vec)))

	dExp_EN<-matrix2dframe(dataExp_EN, nTreatment, sampleSize); #untransformed 
#now do the linear regression with interaction 
	#if(i>=19){
	lreg_EN_lm<-lm(exp~gene*group, data=dExp_EN)#, weight=1/var_pos_vec_all) 
	#lreg_EN<-rlm(exp~gene*group, data=dExp_EN)#, weight=1/var_pos_vec_all) 
	lregSm<-summary(lreg_EN)

	#lreg_EN<-lm(exp~gene*group, data=dExp_EN)#, weight=1/var_pos_vec_all) 
	#lregSm<-summary(lreg_EN)
	#we have the data, what to do.
	#reformat the data from data matrix to 
	#dataframe for the linear regress
	#dExp_EN<-matrix2dframe(dataExp_EN, nTreatment, sampleSize);

	#now do the linear regression with interaction 

	lstToSave<-list("data"=dataExp_EN_List, "lregCoef"=lregSm$coefficients, "var.prior"=r2);
	save(lstToSave, file=paste("result_",i,"_rlm.RData",sep=""));
	#}
}



##first checking for weighted unEqualVar dataset -- Negative control
setwd("E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_IsoConWt")
#-nGene<-2000
#-nTreatment<-2   #number of different beta
#-sampleSize<-5 ##this is the number of repeats for each group
#-alpha.mean<-0  #variance for alpha prior
#-alpha.sigma<-3
#-beta.mean<-0
#-beta.sigma<-2   #variance for beta prior

#-gamma.sigma<-10  # the unscaled factor for gamma given gamma<=>0
#-prob.nonZero<-0.02
#priors for variance distribution
#-d0<-5
#-s0<-2
repeats<-100

#this time of around, we need a fixed variance for the gaussian variance
#epsilon.si<-2 #not use for this case
#call it
set.seed(2005);
#simulate *****the equal variance and !!!!!!!negative control data
#E-equal variance; N-negative control
i<-1
for(i in c(1:repeats))
{
cat("doing ", i, "/",repeats," analyses....\n")
	dataExp_EN_List<-simulateExpression(nGene, nTreatment, sampleSize,
				control.negative=FALSE, control.isotype=TRUE,
				control.index=c(1),
				alpha.mean=alpha.mean, alpha.sigma=alpha.sigma,
				beta.mean=beta.mean, beta.sigma=beta.sigma,
				prob.nonZero=prob.nonZero, gamma.sigma=gamma.sigma,
				#epsilon.si=epsilon.si,
				epsilon.d0=d0, epsilon.s0=s0
				)
#calculate the sample variance from the data
dataExp_EN<-dataExp_EN_List[[1]]
	sMean<-sampleMean(dataExp_EN, nTreatment, sampleSize)
	sVar<-sampleVariance(dataExp_EN, nTreatment, sampleSize)
	sdf<-sampleSizes(dataExp_EN, nTreatment, sampleSize)-1
	#get the sample variance prior by assuming a scale-chisquare
	s2<-as.vector(as.matrix(sVar))
	df<-as.vector(as.matrix(sdf))
	
	r2<-sampleVariancePrior(s2,df)
	s02<-r2[[2]]
	df0<-r2[[1]]
	
	#calculate the posterior variance
	var_pos<-sampleVariancePosterior(sdf,sVar,  df0, s02 )
	var_pos_vec<-(as.vector(t(var_pos)))
	var_pos_vec_all<-rep(var_pos_vec,rep(5, length(var_pos_vec)))

	dExp_EN<-matrix2dframe(dataExp_EN, nTreatment, sampleSize); #untransformed 
#now do the linear regression with interaction 

	lreg_EN<-lm(exp~gene*group, data=dExp_EN, weight=1/var_pos_vec_all) 
	lregSm<-summary(lreg_EN)

	#lreg_EN<-lm(exp~gene*group, data=dExp_EN)#, weight=1/var_pos_vec_all) 
	#lregSm<-summary(lreg_EN)
#we have the data, what to do.
#reformat the data from data matrix to 
#dataframe for the linear regress
#dExp_EN<-matrix2dframe(dataExp_EN, nTreatment, sampleSize);

#now do the linear regression with interaction 

lstToSave<-list("data"=dataExp_EN_List, "lregCoef"=lregSm$coefficients, "var.prior"=r2);
save(lstToSave, file=paste("result_",i,".RData",sep=""));
}


##first checking for weighted unEqualVar dataset -- Negative control
setwd("E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_UnConWt")
#-nGene<-2000
#-nTreatment<-2   #number of different beta
#-sampleSize<-5 ##this is the number of repeats for each group
#-alpha.mean<-0  #variance for alpha prior
#-alpha.sigma<-3
#-beta.mean<-0
#-beta.sigma<-2   #variance for beta prior

#-gamma.sigma<-10  # the unscaled factor for gamma given gamma<=>0
#-prob.nonZero<-0.02
#priors for variance distribution
#-d0<-5
#-s0<-2
repeats<-100

#this time of around, we need a fixed variance for the gaussian variance
#epsilon.si<-2 #not use for this case
#call it
set.seed(2006);
#simulate *****the equal variance and !!!!!!!negative control data
#E-equal variance; N-negative control
i<-1
for(i in c(1:repeats))
{
cat("doing ", i, "/",repeats," analyses....\n")
	dataExp_EN_List<-simulateExpression(nGene, nTreatment, sampleSize,
				control.negative=FALSE, control.isotype=FALSE,
				control.index=c(1),
				alpha.mean=alpha.mean, alpha.sigma=alpha.sigma,
				beta.mean=beta.mean, beta.sigma=beta.sigma,
				prob.nonZero=prob.nonZero, gamma.sigma=gamma.sigma,
				#epsilon.si=epsilon.si,
				epsilon.d0=d0, epsilon.s0=s0
				)
#calculate the sample variance from the data
dataExp_EN<-dataExp_EN_List[[1]]
	sMean<-sampleMean(dataExp_EN, nTreatment, sampleSize)
	sVar<-sampleVariance(dataExp_EN, nTreatment, sampleSize)
	sdf<-sampleSizes(dataExp_EN, nTreatment, sampleSize)-1
	#get the sample variance prior by assuming a scale-chisquare
	s2<-as.vector(as.matrix(sVar))
	df<-as.vector(as.matrix(sdf))
	
	r2<-sampleVariancePrior(s2,df)
	s02<-r2[[2]]
	df0<-r2[[1]]
	
	#calculate the posterior variance
	var_pos<-sampleVariancePosterior(sdf,sVar,  df0, s02 )
	var_pos_vec<-(as.vector(t(var_pos)))
	var_pos_vec_all<-rep(var_pos_vec,rep(5, length(var_pos_vec)))

	dExp_EN<-matrix2dframe(dataExp_EN, nTreatment, sampleSize); #untransformed 
#now do the linear regression with interaction 

	lreg_EN<-lm(exp~gene*group, data=dExp_EN, weight=1/var_pos_vec_all) 
	lregSm<-summary(lreg_EN)

	#lreg_EN<-lm(exp~gene*group, data=dExp_EN)#, weight=1/var_pos_vec_all) 
	#lregSm<-summary(lreg_EN)
#we have the data, what to do.
#reformat the data from data matrix to 
#dataframe for the linear regress
#dExp_EN<-matrix2dframe(dataExp_EN, nTreatment, sampleSize);

#now do the linear regression with interaction 

lstToSave<-list("data"=dataExp_EN_List, "lregCoef"=lregSm$coefficients, "var.prior"=r2);
save(lstToSave, file=paste("result_",i,".RData",sep=""));
}



#================assumption checking======= negative case
#using the simulated data to show the linear regression 
#assumptions before and after the data transformation/variance stabilization
###start doing the simulation
set.seed(2004);
i<-1
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
#now do the transformation
#dataTransformed<-transformData(dataExp_EN,nTreatment, sampleSize)

#with data we can do get population variance for each group
	sMean<-sampleMean(dataExp_EN, nTreatment, sampleSize)
	sVar<-sampleVariance(dataExp_EN, nTreatment, sampleSize)
	sdf<-sampleSizes(dataExp_EN, nTreatment, sampleSize)-1
	#get the sample variance prior by assuming a scale-chisquare
	s2<-as.vector(as.matrix(sVar))
	df<-as.vector(as.matrix(sdf))
	
	r2<-sampleVariancePrior(s2,df)
	s02<-r2[[2]]
	df0<-r2[[1]]
	
	#calculate the posterior variance
	var_pos<-sampleVariancePosterior(sdf,sVar,  df0, s02 )
	var_pos_vec<-(as.vector(t(var_pos)))
	var_pos_vec_all<-rep(var_pos_vec,rep(5, length(var_pos_vec)))

	dExp_EN<-matrix2dframe(dataExp_EN, nTreatment, sampleSize); #untransformed 
#now do the linear regression with interaction 

	wlreg_EN<-lm(exp~gene*group, data=dExp_EN, weight=1/var_pos_vec_all) 
	wlregSm<-summary(wlreg_EN)

	lreg_EN<-lm(exp~gene*group, data=dExp_EN)#, weight=1/var_pos_vec_all) 
	lregSm<-summary(lreg_EN)
###dianostic plots
setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_NegativeCon")
pdf(file="lregDiag.pdf")
op<-par(mfrow=c(2,2))
plot(lreg_EN)
par(op)
dev.off()

#we have the data, what to do.
#reformat the data from data matrix to 
#dataframe for the linear regress
dExp_EN<-matrix2dframe(dataTransformed, nTreatment, sampleSize);#transformed

#now do the linear regression with interaction 

	lreg_EN<-lm(exp~gene*group, data=dExp_EN) 
	lregSm<-summary(lreg_EN)


###dianostic plots
setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_NegativeCon_Trans")
pdf(file="lregDiag.pdf")
op<-par(mfrow=c(2,2))
plot(lreg_EN)
par(op)
dev.off()


#####################done with the assumption checking
##++++++++++++++++++++++++++++++++++++++++

#-------start doing the analysis--------
#setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_NegConWt");

filePath<-"E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_NegConWt";

rtlist<-analyzeData(path=filePath,filename="result_", repeats=100, 
					sample.size=5,
					proportion.nonZero=0.02, 
					object.load="lstToSave",
					mode=1);
					
###adding data together for comparison (see function analyzeData for code)
x_lr_coe_int_q<-cbind(lr_coe_int,lr_coe_int_prob_qvalue)	
xx_lr_coe_int_q<-rbind(c(0,0,0,0,0),x_lr_coe_int_q)
xxx_lr_coe_int_q<-cbind(xx_lr_coe_int_q,prob_ordT_Q)
y_lr_coe_int_q<-cbind(xxx_lr_coe_int_q,roc_response)
 write.table(y_lr_coe_int_q, "./temp.txt", sep="\t")
#####---------------
filePath<-"E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_NegativeCon";

rtlist<-analyzeData(path=filePath,filename="result_", repeats=100, 
					sample.size=5,
					proportion.nonZero=0.02, 
					object.load="lstToSave",
					mode=1);			
	

########=================Isotype control data analysis

#setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_IsotypeCon_Trans");

filePath<-"E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_IsoConWt";

rtlist<-analyzeData(path=filePath,filename="result_", repeats=100, 
					sample.size=5,
					proportion.nonZero=0.02, 
					object.load="lstToSave",
					mode=1);

	
 #######====================No control++++++++++

#setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_UnControl_Trans");

filePath<-"E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_UnControl_Trans";

rtlist<-analyzeData(path=filePath,filename="result_", repeats=100, 
					sample.size=5,
					proportion.nonZero=0.02, 
					object.load="lstToSave_ENC",
					mode=1); # <<----mode of 1
					
					
					
filePath<-"E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_UnControl_Trans";

rtlist<-analyzeData(path=filePath,filename="result_", repeats=100, 
					sample.size=5,
					proportion.nonZero=0.02, 
					object.load="lstToSave_ENC",
					mode=2); #<<---mode of 2
					
	
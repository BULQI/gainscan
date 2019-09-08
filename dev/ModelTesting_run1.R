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
##----------------------------
## ==update===
## Dec 13, 2016 Feng
##   add code to do FDR so to see how it behave with FDR analysis
##			reasoning: why do this? Since don't know the prob of differential interaction
##
##  =======update========
## Dec 15, 2016 Feng
##    write the code to wrap up analysis into functions
##		see functionAnalysis.r for the function definition
##		
##
##########################################
library(ARPPA)
library(qvalue)
#install.packages("pROC")
library(pROC)
##first checking for homoscedastical dataset
setwd("E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\EqualVar_NegativeCon")
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


#this time of around, we need a fixed variance for the gaussian variance
epsilon.si<-2
#call it
set.seed(2004);
#simulate *****the equal variance and !!!!!!!negative control data
#E-equal variance; N-negative control
repeats<-100
for(i in c(1:repeats))
{
cat("doing ", i, "/",repeats," analyses....\n")
	dataExp_EN_List<-simulateExpression(nGene, nTreatment, sampleSize,
				control.negative=TRUE, control.isotype=FALSE,
				control.index=c(1),
				alpha.mean=alpha.mean, alpha.sigma=alpha.sigma,
				beta.mean=beta.mean, beta.sigma=beta.sigma,
				prob.nonZero=prob.nonZero, gamma.sigma=gamma.sigma,
				epsilon.si=epsilon.si #epsilon.d0=d0, epsilon.s0=s0
				)
#calculate the sample variance from the data
dataExp_EN<-dataExp_EN_List[[1]]
svList<-sampleVariance(dataExp_EN,nTreatment, sampleSize)
#reformat it into a vector
s2<-as.vector(as.matrix(svList[[1]]))
df<-as.vector(as.matrix(svList[[2]]))
r2<-sampleVariancePrior(s2,df)

#we have the data, what to do.
#reformat the data from data matrix to 
#dataframe for the linear regress
dExp_EN<-matrix2dframe(dataExp_EN, nTreatment, sampleSize);

#now do the linear regression with interaction 
lreg_EN<-lm(exp~gene*group, data=dExp_EN)
lregSm<-summary(lreg_EN)
lstToSave<-list("data"=dataExp_EN_List, "lregCoef"=lregSm$coefficients, "var.prior"=r2);
save(lstToSave, file=paste("result_",i,".RData",sep=""));
}

###############
set.seed(2005)
setwd("E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\EqualVar_IsotypeCon")
for(j in c(2:repeats))
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
				epsilon.si=epsilon.si #epsilon.d0=d0, epsilon.s0=s0
				)
dataExp_EI<-dataExp_EI_List[[1]]
#calculate the sample variance from the data
svList<-sampleVariance(dataExp_EI,nTreatment, sampleSize)
#reformat it into a vector
s2<-as.vector(as.matrix(svList[[1]]))
df<-as.vector(as.matrix(svList[[2]]))
r2<-sampleVariancePrior(s2,df)

#we have the data, what to do.
#reformat the data from data matrix to 
#dataframe for the linear regress
dExp_EI<-matrix2dframe(dataExp_EI, nTreatment, sampleSize);

#now do the linear regression with interaction 
lreg_EI<-lm(exp~gene*group, data=dExp_EI)
lregSm_EI<-summary(lreg_EI)
lstToSave_EI<-list("data"=dataExp_EI_List, "lregCoef"=lregSm_EI$coefficients, "var.prior"=r2);
save(lstToSave_EI, file=paste("result_",j,".RData",sep=""));
cat("Done!!\n")
}

#####################
set.seed(2006)
setwd("E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\EqualVar_Uncontrol")
for(k in c(2:repeats))
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
				epsilon.si=epsilon.si #epsilon.d0=d0, epsilon.s0=s0
				)
dataExp_ENC<-dataExp_ENC_List[[1]]
#calculate the sample variance from the data
svList<-sampleVariance(dataExp_ENC,nTreatment, sampleSize)
#reformat it into a vector
s2<-as.vector(as.matrix(svList[[1]]))
df<-as.vector(as.matrix(svList[[2]]))
r2<-sampleVariancePrior(s2,df)

#we have the data, what to do.
#reformat the data from data matrix to 
#dataframe for the linear regress
dExp_ENC<-matrix2dframe(dataExp_ENC, nTreatment, sampleSize);

#now do the linear regression with interaction 
lreg_ENC<-lm(exp~gene*group, data=dExp_ENC)
lregSm_ENC<-summary(lreg_ENC)
lstToSave_ENC<-list("data"=dataExp_ENC_List, "lregCoef"=lregSm_ENC$coefficients, "var.prior"=r2);
save(lstToSave_ENC, file=paste("result_",k,".RData",sep=""));
cat("Done!!\n")
}
#now summary the data first
aggregate(exp~gene*group, dExp_ENC, FUN=mean)

###dianostic plots
op<-par(mfrow=c(2,2))
plot(lreg_EN)
par(op)

#================analysis= negative control======
###start doing the summary statistics of the repeated runs
#
#loading the data
repeats<-10
#setwd("E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\EqualVar_NegativeCon");
filePath<-"E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\EqualVar_NegativeCon"

rtlist<-analyzeData(path=filePath,filename="result_", repeats=100, 
					sample.size=5,
					proportion.nonZero=0.02);
					

########=================Isotype control data analysis
#loading the data

setwd("E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\EqualVar_IsotypeCon");
#setwd("~/Desktop/arppr/EqualVar_NegativeCon/")
filePath<-"E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\EqualVar_IsotypeCon"

rtlist<-analyzeData(path=filePath,filename="result_", repeats=100, 
					sample.size=5,
					proportion.nonZero=0.02,
					object.load="lstToSave_EI",
					mode=1);
 
 #######====================No control++++++++++

#loading the data
repeats<-100
#setwd("E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\EqualVar_UnControl");
#setwd("~/Desktop/arppr/EqualVar_NegativeCon/")

filePath<-"E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\EqualVar_UnControl"

rtlist<-analyzeData(path=filePath,filename="result_", repeats=100, 
					sample.size=5,
					proportion.nonZero=0.02,
					object.load="lstToSave_ENC", 
					mode=1);
					
rtlist<-analyzeData(path=filePath,filename="result_", repeats=100, 
					sample.size=5,
					proportion.nonZero=0.02,
					object.load="lstToSave_ENC", 
					mode=2);


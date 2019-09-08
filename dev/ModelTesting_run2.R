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
######
##  update, Dec 16, 2016 feng@bu
##	change the code of analysis to functions defined @ functionAnlysis.r
##   for the old code, see ModelTesting_run2.R.old
##########################################
library(ARPPA)
##first checking for homoscedastical dataset
setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_NegativeCon")
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

#we have the data, what to do.
#reformat the data from data matrix to 
#dataframe for the linear regress
dExp_EN<-matrix2dframe(dataExp_EN, nTreatment, sampleSize);

#now do the linear regression with interaction 
#now try to read the gpr files
 gpr <- "H:\\feng\\LAB\\MSI\\AbSynthesis\\proteinArray\\ProtoArray_20161026_isotype"
 targets <- list.files(path="H:\\feng\\LAB\\hg\\proteinArray_Masa\\U19IOF",
					pattern = "^Isotype1026", full.names = TRUE)
 ad.elist <- loadGPR(gpr.path=gpr, targets.path=targets, array.type="ProtoArray",aggregation="none")

#save data
 save(ad.elist, file=paste(gpr, "/AD.RData",sep=""), compress="xz")
#load data
 load(paste(gpr,"/AD.RData",sep=""))
 
 
 #now plot the raw data
 plot(density(ad.elist$E[,1]), main="raw data distribution", xlab="intensity", ylab="freq", col="red", type="p")
 lines(density(ad.elist$E[,2]), col="green")
 
 plot(density(ad.elist$C[,1]))
lregSm<-summary(lreg_EN)
lstToSave<-list("data"=dataExp_EN_List, "lregCoef"=lregSm$coefficients, "var.prior"=r2);
save(lstToSave, file=paste("result_",i,".RData",sep=""));
}

###############
set.seed(2005)
setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_IsotypeCon")
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
setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_UnControl")
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

#================data analysis=======Negative control
###start doing the summary statistics of the repeated runs

#setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_NegativeCon");

filePath<-"E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_NegativeCon"

rtlist<-analyzeData(path=filePath,filename="result_", repeats=100, 
					sample.size=5,
					proportion.nonZero=0.02, 
					object.load="lstToSave",
					mode=1);


########=================Isotype control data analysis
#loading the data
#setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_IsotypeCon");

filePath<-"E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_IsotypeCon";

rtlist<-analyzeData(path=filePath,filename="result_", repeats=100, 
					sample.size=5,
					proportion.nonZero=0.02, 
					object.load="lstToSave_EI",
					mode=1);
					

 
 #######====================No control++++++++++

#setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_UnControl");

filePath<-"E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_UnControl";

rtlist<-analyzeData(path=filePath,filename="result_", repeats=100, 
					sample.size=5,
					proportion.nonZero=0.02, 
					object.load="lstToSave_ENC",
					mode=1);
					

filePath<-"E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_UnControl";

rtlist<-analyzeData(path=filePath,filename="result_", repeats=100, 
					sample.size=5,
					proportion.nonZero=0.02, 
					object.load="lstToSave_ENC",
					mode=2);



##R code
##Tutorial_v1.2.R
##this R code is on top of Tutorial_v1.0
##in here we run fitting of arrays profiling autoreactive antibodies
##by Feng @ 1/28/2019
##
#------updated to do 
#==no BC data analysis and fitting
#---3/31/2019
#-=======++++++
##====================

#load libraries
library(PAA)
library(PAST)
library(limma)
library(minpack.lm)

#now start reading the data 
 pmt_saturated<-2^16-1
gpr <- "E:\\feng\\LAB\\MSI\\AbSynthesis\\proteinArray\\protoArray\\11_28_2016\\11_28_2016_gpr"
 targets <- list.files(path="E:\\feng\\LAB\\MSI\\AbSynthesis\\proteinArray\\protoArray\\11_28_2016\\11_28_2016_gpr",
					pattern = "^target_6192_1182.txt", full.names = TRUE)
#read the saved data
 load(paste(gpr,"/AD_6192_1182.RData",sep=""))
 pmts.two<-c(475,525,575,600)

 #try to load .
 setwd( gpr)
 
 load("AD_BC_6192_1182.RData")
x<-pmts.two 
  
  #+++++++++now we do no BC fitting +++++++++
  ad.elistBC.two<-ad.elist.two
  #++++++++++++++++++++++++++++++++++
  
  
#get the testing data by block for testing 
#now getting data ready to do the fitting......
		object<-ad.elistBC.two;
#object<-ad.elistBC; x<-x;G<-pmt_saturated; block.size<-4; fit.mode<-"control"; residual.mode="log";debug=T
p_elists<-gainAdjust_fit_Pbp(object=ad.elistBC.two, x=x, 
								#ylong=F, aggregated=F,
								 B=log(30), 
								G=pmt_saturated, #phi, delta,
								block.size=4, fit.mode="control",
								 residual.mode="log",
								debug=T )
p_elists.auto<-p_elists
								
setwd("e:/feng/LAB/mammr/2018Jan_ProgressReport/ProtoArray/Figure5_SignalToNoiseRatio/")
#save(p_elists.auto, file="PMT_Pbp_p_elist_autoreactive.RData")
save(p_elists.auto, file="PMT_Pbp_p_elist_autoreactive_noBC.RData") #<----to do no BC and with log(x) residual 3/31


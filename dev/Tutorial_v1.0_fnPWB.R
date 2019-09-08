##R code
##Tutorial_v1.0.R
##this R code is used to test and try the new version of PAST, 
##which fits the data with the PMT gain power function with baseline.
##version 1.0
##
##by Feng @ 1/28/2019
##
#++++++++++++updated 
# 3/31/2019
#++to do no BC data analysis and fitting
#_______________

##====================

#load libraries
library(PAA)
library(PAST)
library(limma)
library(minpack.lm)

#now start reading the data 
 pmt_saturated<-2^16-1
 #now try to read the gpr files
 gpr <- "E:\\feng\\LAB\\MSI\\AbSynthesis\\proteinArray\\protoArray\\11_28_2016\\iggisotyprctrl_feng_12052016"
  #gpr <- "/home/feng/Windows/D/feng/LAB/MSI/AbSynthesis/proteinArray/protoArray/11_28_2016/iggisotyprctrl_feng_12052016"
 targets <- list.files(path=gpr,
					pattern = "^target_second.txt", full.names = TRUE)
#ad.elist <- loadGPR(gpr.path=gpr, targets.path=targets, array.type="ProtoArray",aggregation="none")
		#save(ad.elistBC,array1.C,array1.E, file= "AD_BC.RData", compress="xz")
#load data
 load(paste(gpr,"/AD_Acquire5th.RData",sep=""))
 
 ######==========start doing the 5PL fitting 
 ###visualize the data first to see whether it is a good model
 pmts<-seq(250, 650, by=50)
 
 #background correction
 #setwd( "E:\\feng\\LAB\\MSI\\AbSynthesis\\proteinArray\\protoArray\\11_28_2016\\iggisotyprctrl_feng_12052016")
 setwd(gpr);
 #ad.elistBC <- bc(ad.elist, method="normexp",
 # normexp.method="saddle")
 
	#save(ad.elistBC, file= "AD_BC_6thPAST.RData", compress="xz")
 #load data
 load("AD_BC_6thPAST.RData")
  x<-pmts 

aggregated<-F


#....<----+++++++++++++++++++++++++
#+++++++++now we fit with no BC data+++++++++
ad.elistBC<-ad.elist;
#++++++++++++++++++++++++++++++++++++



#get the testing data by block for testing 
#now getting data ready to do the fitting......
		object<-ad.elistBC
#now we have the data read in. 
#for the testing at this stage, let's process only one array.
index.firstArr<-c(1:9)
object$E<-object$E[,index.firstArr];
object$C<-object$C[,index.firstArr];
object$targets<-object$targets[index.firstArr,];
object.oneArray<-object
#now getting a blocked data to do fitting.
num.block<-24
block.start<-1
cgenes<-object$cgenes;
genes<-object$genes;

index.block<-which(is.element(cgenes[,"Block"],c(block.start:(block.start+num.block-1))))
object$C<-object$C[index.block,]	
object$cgenes<-object$cgenes[index.block,]

index.block<-which(is.element(genes[,"Block"],c(block.start:(block.start+num.block-1))))
object$E<-object$E[index.block,]
object$genes<-object$genes[index.block,]


#now let's feed in to test
inits<-gainAdjust_prepareInits_Pbp(y=object$C, x, aggregated=F, data.log=F)

#now first testing the power function with baseline.
#B<-20;
#Ys<-f_PowerBaselinePMT(c(B, 2^16-1, inits$inits.phi[1], inits$inits.delta),x);

B<-log(1.1)
pars<-c(B, inits$inits.delta, inits$inits.phi)
y<-object$C
y<-gainAdjust.aggregate(y, mode="geometric")
G<-2^16-1
residual.mode<-"log"

#rs<-gainAdjust_fnResPbp(pars, y, x, residual.mode="log")

#ft<-nls.lm(par=pars, fn=gainAdjust_fnResPbp,
#					#y=log(yinput.aligned), 
#					y=y, 
#					x=x, 
#					residual.mode="log",
#					control = nls.lm.control(nprint=1, maxit=100)
#				)

#save(file="ft.Rdata", ft);
#c(max(ft$par[-c(1,2)]),min(ft$par[-c(1,2)]))
#sink(file="temp.txt")
system.time(
ftf<-gainAdjust_fit_Pbp_data(object=object$C, x=x, 
								#ylong=F,
								data.aggregated=F,
								 B=B, G=2^16-1, #phi=log(inits$inits.phi), delta=inits$inits.delta,
								 residual.mode="log",
								debug=T )
							)
c(max(ftf$par[-c(1,2)]),min(ftf$par[-c(1,2)]))							
#sink();
	
#now plot the fitting.
gainAdjust_plot_Pbp_raw(LMfit=ftf,y, x, data.aggregated=T,
			#model="fix.both",
			#B,
			G=pmt_saturated, log="xy",
			#delta, phi, ref.index=1,
			filename=NULL);
			
#now plot them
#y=ydata.blocked; x=x;LMfit=ftf;data.aggregated=T;G=pmt_saturated;ref.index=1;log="xy";type="l"
gainAdjust_plot_Pbp(LMfit=ftf,y, x, data.aggregated=T,
			G=pmt_saturated, 
			ref.index=1,  #used to indicate which line is used as the reference, all other ones is shifted towards it.
			log="xy", type="l",
			filename=NULL);

#now testing the fit the phi with fixed B and delta

y_E<-gainAdjust.aggregate(object$E, mode="geometric")
inits<-gainAdjust_prepareInits_Pbp(y=y_E, x, aggregated=T, data.log=F);
#test the first series
rs<-gainAdjust_fnResPbp_phi(pars=log(inits$inits.phi[1]), 
		y=y_E[1,], #the response, assuming log transformed
		x, #the independent variable, log-transformed
		B=ftf$par[1], G=pmt_saturated, delta=ftf$par[2],
		residual.mode="log"
		#xmid, #x value at the middel point of y
		#scal, #slope at xmid
		#g #the asymetric factor
		)

#object=y_E;x=x; B=ftf$par[1]; G=pmt_saturated; delta=ftf$par[2];residual.mode="log"; debug=T;
f_E<-gainAdjust_fit_Pbp_phi(object=y_E,x=x, B=ftf$par[1], G=pmt_saturated, delta=ftf$par[2],residual.mode="log", debug=T);

#now testing the fitting for one array by blocks.
#G<-pmt_saturated;block.size<-4;fit.mode<-"control";residual.mode<-"log";debug<-T;object<-object.oneArray
p_elist<-gainAdjust_fit_Pbp_array(object=object.oneArray, x=x, 
								G=pmt_saturated, B=log(30),
								block.size=4, fit.mode="control",
								 residual.mode="log",
								debug=T )
#object<-ad.elistBC; x<-x;G<-pmt_saturated; block.size<-4; fit.mode<-"control"; residual.mode="log";debug=T
p_elists<-gainAdjust_fit_Pbp(object=ad.elistBC, x=x, 
								#ylong=F, aggregated=F,
								# B=pmt_baseline, 
								B=log(30),
								G=pmt_saturated, #phi, delta,
								block.size=4, fit.mode="control",
								 residual.mode="log",
								debug=T )
								
setwd("e:/feng/LAB/mammr/2018Jan_ProgressReport/ProtoArray/figure3_variance/")
#save(p_elists, file="PMT_Pbp_p_elist_Isotype_lowerLimit.RData")
save(p_elists, file="PMT_Pbp_p_elist_Isotype_lowerLimit_noBC.RData") #,-this is updated 3/31 using no BC and fixed other things, with
				#log(x) scaled residual 


		#		++++++++++++---------------
		#-------------------
#now testing the ill cases. almost 
fix_B=  log(2.189615) 

fix_delta=6.877512
input_x<-c(700,550)
input_y<-c( 45.34970,             9.238788,15.27166,            11.713661)

ft<-nls.lm(log(0.005), fn=gainAdjust_fnResPbp_phi, #lower=lower, 
						y=input_y, x=input_x, B=fix_B, G=G, delta=fix_delta, residual.mode="log",
						control = nls.lm.control(nprint=1, maxit=100)
						);
						
#testing the initials
#y=matrix(input_y, nrow=2,byrow=T); x=input_x; aggregated=T; data.log=F;
inits<-gainAdjust_prepareInits_Pbp(y=matrix(input_y, nrow=2, byrow=T), x=input_x, aggregated=T, data.log=F);

#############the following part is experimenting the fitting of the Anti-histone antibody.
#read data
setwd("E:\\feng\\LAB\\MSI\\AbSynthesis\\proteinArray\\20180911_exp_results");
#setwd("C:\\feng\\LAB\\MSI\\AbSynthesis\\proteinArray\\20180911_exp_results");

pmt_saturated<-2^16-1

#first you need to create a target file including all the meta data for each array, like below
targets <- read.table(file=list.files("./",
					pattern = "^target.txt", full.names = TRUE), header=TRUE)
# print(targets[1:2,])
 
 #now try to read the gpr files
 gpr <- "E:\\feng\\LAB\\MSI\\AbSynthesis\\proteinArray\\20180911_exp_results"
 #gpr <- "C:\\feng\\LAB\\MSI\\AbSynthesis\\proteinArray\\20180911_exp_results"
 targets<-list.files("./",
					pattern = "^target.txt", full.names = TRUE)
 
 load(paste(gpr,"/AD.RData",sep=""))
 
 pmts<-c(350,450,550,600,650,700);
load(paste(gpr,"/AD_BC.RData",sep=""));
	ad.elistBC<-ad.elist #<---------now doing raw no BC'ed data
	ad.elistBC.hf<-ad.elistBC 
	ad.elistBC.hf$E<-ad.elistBC$E[,c(19:36,1:18)]
	ad.elistBC.hf$C<-ad.elistBC$C[,c(19:36,1:18)]
	ad.elistBC.hf$targets<-ad.elistBC$targets[c(19:36,1:18),]

	x<-pmts
	object<-ad.elistBC.hf
#now get one array to do testing.
index.firstArr<-c(1:6)
object$E<-object$E[,index.firstArr];
object$C<-object$C[,index.firstArr];
object$targets<-object$targets[index.firstArr,];
object.oneArray<-object
#now getting a blocked data to do fitting.
num.block<-4
block.start<-9
cgenes<-object$cgenes;
genes<-object$genes;

index.block<-which(is.element(cgenes[,"Block"],c(block.start:(block.start+num.block-1))))
object$C<-object$C[index.block,]	
object$cgenes<-object$cgenes[index.block,]

index.block<-which(is.element(genes[,"Block"],c(block.start:(block.start+num.block-1))))
object$E<-object$E[index.block,]
object$genes<-object$genes[index.block,]

B<-log(30)
G<-2^16-1
system.time(
#object=object$C; x=x; B=B; G=2^16-1; residual.mode="log";debug=T;data.aggregated=F 
ftf<-gainAdjust_fit_Pbp_data(object=object$C, x=x, 
								#ylong=F,
								 data.aggregated=F,
								 B=B, G=2^16-1, #phi=log(inits$inits.phi), delta=inits$inits.delta,
								 residual.mode="log",#minimum.phi=-7,
								debug=T )
							)
#y=y; x; aggregated=TRUE; data.log=F

gainAdjust_fit_Pbp_phi(object=y[62,],x=x, B=ft$par[1], delta=ft$par[2])
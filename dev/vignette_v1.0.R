##R code
##v1.0 for vignette
##this R code is used to test and imbeded into the vignette of gainscan, 
##which fits the data with the PMT gain power function with baseline.
##version 1.0
##
##by Feng @ 9/28/2019
##
#++++++++++++updated 
# originally copied from the tutorial v1.0
#_______________

##====================

#load libraries

library(gainscan)
library(limma)
library(minpack.lm)

#now start reading the data 
 pmt_saturated<-2^16-1
 #now try to read the gpr files
 gpr<-system.file("extdata", package="gainscan")
 targets <- list.files(path=gpr,
					pattern = "target.txt", full.names = TRUE)
tbl.targets<-read.table(targets, header=T,sep="\t")
print(tbl.targets)

#->don't use this this is from PAA (use the one from the gainscan package) ad.elist <- loadGPR(gpr.path=gpr, targets.path=targets, array.type="ProtoArray",aggregation="none")
setwd(gpr);
adata.elist <- importGPR(dir.GPR=gpr, design.file=targets, type="ProtoArray")

#save(ad.elistBC,array1.C,array1.E, file= "AD_BC.RData", compress="xz")
#load data
# load(paste(gpr,"/AD_Acquire5th.RData",sep=""))
 
 ######==========start doing the 5PL fitting 
 ###visualize the data first to see whether it is a good model
 pmts<-seq(250, 650, by=50)
 
 #background correction
 #setwd( "E:\\feng\\LAB\\MSI\\AbSynthesis\\proteinArray\\protoArray\\11_28_2016\\iggisotyprctrl_feng_12052016")
 #setwd(gpr);
 adata.elistBC <- bc(adata.elist, method="normexp",
  normexp.method="saddle")

  #compare expression data
 cbind("before"=adata.elist$E[1:8,1], "after"=adata.elistBC$E[1:8,1])
	#save(ad.elistBC, file= "AD_BC_6thPAST.RData", compress="xz")
 #compare control data
  cbind("before"=adata.elist$C[1:8,1], "after"=adata.elistBC$C[1:8,1])
 

B<-log(min(adata.elist$C)); #this has to be as close as possible, otherwise might affect the fitting.

#object<-ad.elistBC; x<-x;G<-pmt_saturated; block.size<-4; fit.mode<-"control"; residual.mode="log";debug=T
p_elists<-gainAdjust_fit_Pbp(object=adata.elist, x=pmts, 
								#ylong=F, aggregated=F,
								# B=pmt_baseline, 
								B=B,
								G=pmt_saturated, #phi, delta,
								block.size=2, fit.mode="control",
								 residual.mode="log",
								debug=F )
								
#doing plot of the density
op<-par(
			mfrow=c(2,1),
			pty="m"
		)
 
 plot(density(log(adata.elist$E[,17])),main="Distribution of array signal intensities (raw)",
	xlab="intensity (log)",
	)
 lines(density(log(adata.elist$E[,8])))
 lines(density(log(adata.elist$E[,26])))
 
 plot(density(p_elists$elist$E[,1]), main="Distribution of array signal intensities (fitted)",
	xlab="incide light signal (log)",
	)
 lines(density(p_elists$elist$E[,2]))
 lines(density(p_elists$elist$E[,3]))
 
par(op)

#start plotting the distribution of signal intensity
object.pbp<-scaleByDelta(p_elists$elist, mean(p_elists$delta))  #<-scale it, not necessary, but just try
#now do the RLM normalization
object.pbp.norm<-normalizeArrays.RLM(data=object.pbp,controls="Anti-HumanIgG",
				method="RLM", #only implemented RLM for now
				log=F, log.base=exp(1), coding="Deviation"
			)

	object.raw.600<-adata.elist
	object.raw.600$E<-adata.elist$E[,c(8,17,26)];
	object.raw.600$C<-adata.elist$C[,c(8,17,26)];
	object.raw.norm<-normalizeArrays.RLM(data=object.raw.600,controls="Anti-HumanIgG",
				method="RLM", #only implemented RLM for now
				log=TRUE, log.base=exp(1), coding="Deviation"
			)
#doing plot of the density
op<-par(
			mfrow=c(2,1),
				pty="m",
				mar=c(4,4,0.1,0.5),
				mgp=c(1.5,0.5,0)
		)
 
 plot(density((object.raw.norm$E[,3])),main="",
	xlab="intensity (log)",
	)
 lines(density((object.raw.norm$E[,1])))
 lines(density((object.raw.norm$E[,2])))
 legend(x=5.5,y=1.0, legend="signal intensities (raw+RLM)",lty=1)
 
 plot(density(object.pbp.norm$E[,1]), main="",
	xlab="incide light signal (log)",
	)
 lines(density(object.pbp.norm$E[,2]))
 lines(density(object.pbp.norm$E[,3]))
 legend(x=-38.5,y=0.6, legend="after intensities (gainScan+RLM)",lty=1)
par(op)

#calculate the interarray CVs
	#for raw data with PMT600
	dtf.norm<-array.aggregate(object.raw.norm, log=F, method="Arithmetic")
	
	cDtf.norm<-dtf.norm$C
	cgDtf.norm<-dtf.norm$cgenes
	interArray_ccv.norm<-interArrayCVs(cDtf.norm, cgDtf.norm)#sqrt(apply(cDtf.norm, 1, var))/apply(cDtf.norm, 1, 
	
	#now do the PMT fitted and normalizatio 
	cDtf.pbp.norm<-object.pbp.norm$C
	cgDtf.pbp.norm<-object.pbp.norm$cgenes
	interArray_ccv.pbp.norm<-interArrayCVs(cDtf.pbp.norm, cgDtf.pbp.norm)#sqrt(apply(cDtf.5pl.norm, 1,		
	
	#show the mean CVs
	mean(interArray_ccv.norm$cv)#0.052
	mean(interArray_ccv.pbp.norm$cv)

	#start plotting boxplot 					
	interArray_ccv.norm$type<-"RLM Norm"
	ccvByType<-interArray_ccv.norm
	
	interArray_ccv.pbp.norm$type<-"Gain-Scan+RLM"
	ccvByType<-rbind(ccvByType,interArray_ccv.pbp.norm)
	
	ccvByType$type<-factor(ccvByType$type,levels<-c("RLM Norm", "Gain-Scan+RLM"))
	
	boxplot(ccvByType[,1]~type, data=ccvByType,log="",
		main=("InterArray variance"),
		xlab="",ylab="CV", 
		ylim=c(0,0.065),#max(ccvByType[,1])
		#),
		col=2, outline=FALSE
		);

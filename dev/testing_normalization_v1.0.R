#code to do analysis of first run of protoArray with isotype control
#for calculating
#Update
#------------
#12/11/17
#fit the 5PL using sboner. Not sure whether it is going to work. Just try.
##############

# only if you install a Bioconductor package for the first time
# source("http://www.bioconductor.org/biocLite.R")

#load library
# library("BiocInstaller")
# biocLite("PAA", dependencies=TRUE) #install
#
#load the library

library(PAA)
library(limma)
library(minpack.lm)
library(PAST)
library(ggplot2) 

	#now try to load the dataset from Sboner 2009
	#now try to read the gpr files
	gpr_sboner <- "E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\sbnorRLM_data\\ProteinChipRLM_data"
	#targets_sboner <- list.files(path="E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\sbnorRLM_data\\ProteinChipRLM_data",
	#				pattern = "^target.txt", full.names = TRUE) ##this one is to read all the arrays
	targets_sboner <- list.files(path="E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\sbnorRLM_data\\ProteinChipRLM_data",
					pattern = "^target_PCSD1.txt", full.names = TRUE) # to read only PCS 100 and D1 100
	#sboner.elist <- loadGPR(gpr.path=gpr_sboner, targets.path=targets_sboner,
	#	array.type="ProtoArray",aggregation="none",
	#	description="Name", description.features="^Hs~", description.discard="Empty")

	setwd("E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\sbnorRLM_data\\ProteinChipRLM_data");
	#save data
	#save(sboner.elist, file=paste(gpr_sboner, "/sboner_PCSD1.RData",sep=""), compress="xz")
	#load data
	#load(paste(gpr_sboner,"/sboner_PCSD1.RData",sep=""))
 
	###now add background correction
	#sboner.elistBC <- bc(sboner.elist, method="normexp", normexp.method="saddle")
	#save(sboner.elistBC, file=paste(gpr_sboner,"/sboner_PCSD1_BC.RData", sep=""), compress="xz")
	load(paste(gpr_sboner,"/sboner_PCSD1_BC.RData", sep=""))
	#now we get data ready for 5pl fitting
	#we only need one array for this testing.
	
	sboner_HiLow_TwoArray<-sboner.elistBC
	sboner_HiLow_TwoArray$C<-sboner.elistBC$C[,c(1,3)]
	sboner_HiLow_TwoArray$E<-sboner.elistBC$E[,c(1,3)]
	sboner_HiLow_TwoArray$targets<-sboner.elistBC$targets[c(1,3),]	
	
	#start doing the pmt 5pl fitting
	data<-sboner_HiLow_TwoArray
	
	#now testing data frame reformatting
	cdtf<-cmatrix2dataframe(data$C, data$cgenes)
	cdtf.blank<-cmatrix2dataframe(data$C,data$cgenes,control="Blank", log=TRUE)
	cdtf.IgG1<-cmatrix2dataframe(data$C,data$cgenes,control="HumanIgG1", log=TRUE)
	cdtf.IgG<-cmatrix2dataframe(data$C,data$cgenes,control="HumanIgG", log=TRUE)
	
	fit.RLM.C<-fit.RLM(cdtf.IgG);
	
	#start testing read the fitting parameters
	factors<-c("array","block");
	factor.levels<-c(2,48);
	factor.first<-c(colnames(data$C)[1],unique(data$cgenes$Block)[1])
	fit<-fit.RLM.C
	fct.effects<-readFactors(fit,factors, factor.levels, factor.first);
	
	#testing normalization
	control<-"HumanIgG";
	log<-TRUE;
	log.base<-exp(1);
	
	norm.elist<-normalizeArrays.RLM(data=data,controls="HumanIg",
				method="RLM", #only implemented RLM for now
				log=TRUE, log.base=2
			)
	#now for checking the results. first the unchanged
	
	b<-3
	for(b in 1:length(fct.effects[["block"]]))
	{
	cat(paste("doing the ",b, "th block.......",sep=""));
	index.block<-which(data$genes$Block==b)
	
	
	#check for block effects
	effs.obs<-log(data$E[index.block,])-norm.elist$E[index.block,]-matrix(rep(fct.effects[["array"]],length(index.block)),nrow=length(index.block),ncol=length(fct.effects[["array"]]),byrow=T)
	effs.expected<- matrix(fct.effects[["block"]][b],nrow=length(index.block),ncol=dim(data$E)[2],byrow=T)
	
	epsilon<- 1e-12
	stopifnot(abs(effs.obs-effs.expected)<1e-12)
	cat("Done!!!!!!!!!!\n");
	}
	#check for array effects.
	a<-2
	index.array<-a
	effs.obs<-log(data$E[,index.array])-norm.elist$E[,index.array]-fct.effects[["block"]][as.character(data$genes$Block)]
	effs.expected<-matrix(rep(fct.effects[["array"]][index.array],dim(data$E)[1]),nrow=dim(data$E)[1],ncol=length(index.array),byrow=T)
	
	stopifnot(abs(effs.obs-effs.expected)<1e-12)
	
	####now compare with PAA normalization
	###now normalization using PAA package
	sboner_norm<-normalizeArrays(elist=data, method="rlm",
		cyclicloess.method="fast", control ="HumanIgG")
	
	####array effects
	norm.elist$E[11100:11110,1]-norm.elist$E[11100:11110,2]
	#[1] 0.1515504 1.8480025
	sboner_norm$E[11100:11110,1]-sboner_norm$E[11100:11110,2]
	#[1] 0.1511731 1.8476252

	##block effects
	norm.elist$E[1:10,1]-norm.elist$E[11385:11394,1]
	sboner_norm$E[1:10,1]-sboner_norm$E[11385:11394,1]
	
	#now testing the simple coding
	coding<-"Simple"
	data<-cdtf.IgG
	
	#f.simple<-fit.RLM(cdtf.IgG,coding="Simple")
	f.dummy.old<-fit.RLM(cdtf.IgG,coding="Dummy"); 
	ndata.dummy<-normalizeArrays.RLM(data=data,controls="HumanIg",
				method="RLM", #only implemented RLM for now
				log=TRUE, log.base=2, coding="Dummy"
			)
	ndata.simple<-normalizeArrays.RLM(data=data,controls="HumanIg",
				method="RLM", #only implemented RLM for now
				log=TRUE, log.base=2, coding="Simple"
			)
	#this is is old one, hasn't been fixed with new code to do normalization based on coding
	i<-1
	ndata.simple$E[i:(i+20),1]-ndata.simple$E[i:(i+20),2]
	ndata.dummy$E[i:(i+20),1]-ndata.dummy$E[i:(i+20),2]
	#[1] 0.1515504 1.8480025
	sboner_norm$E[i:(i+20),1]-sboner_norm$E[i:(i+20),2]
	#[1] 0.1511731 1.8476252
	
	#----control
	data_c<-data
	data_c$E<-data$C
	data_c$cgene<-data$cgenes
	sboner_norm_c<-normalizeArrays(elist=data_c, method="rlm",
		cyclicloess.method="fast", control ="HumanIgG")
	
	ndata.simple$C[i:(i+20),1]-ndata.simple$C[i:(i+20),2]
	ndata.dummy$C[i:(i+20),1]-ndata.dummy$C[i:(i+20),2]
	#[1] 0.1515504 1.8480025
	sboner_norm_c$E[i:(i+20),1]-sboner_norm_c$E[i:(i+20),2]
	#[1] 0.1511731 1.8476252

	
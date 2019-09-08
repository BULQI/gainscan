###Analysis to see the model assumptions.
# Question to ask ? is it necessary to run a complex model or a simple linear model??
# 
#   Y=a+b+e or Y=a+b+c*e   is c necessary? is this a constant error or not??
#
# ANSWER: 1.definitely the array effects due to primary or secondary antibodies
#		are multiplicative instead of simply additive. Please see the results
#		in Analysis_proteinWithControl.r 
# 		2. the block effects are checked in this analysis.
#
# only if you install a Bioconductor package for the first time
# source("http://www.bioconductor.org/biocLite.R")
# # else
# library("BiocInstaller")
# biocLite("PAA", dependencies=TRUE)
 
 #library(PAA)

library(ARPPA) #this is the my own package to run analysis

#now start reading the text exported data
#this is the batch done on 10/08/2015
datapath<-system.file("extdata", package="ARPPA")
targets <- list.files(system.file("extdata", package="ARPPA"),
	 pattern = "targets_text_Batch1", full.names=TRUE) 
elist2<-importTextData(dataFilePath=datapath, targetFile=targets, start.data=51, nrows.data=18803-53,
			start.control=18803, nrows.control=23286,aggregation="geoMean",
			as.is=TRUE, header=TRUE,sep="\t",na.strings="", quote="")

#===>background correction
library(limma)

elist2 <- backgroundCorrect(elist2, method="normexp",
 normexp.method="saddle")
 
 #now we need to rearrange the control and testing proteins to background correct control
 #it works this way that background correction only working on elist$E
 #so we have to do the rearrangement
 elist_c<-elist2
 elist_c$E<-elist2$C
 elist_c$Eb<-elist2$Cb
 #this is good enough, we don't have to change elistC$C, since other fields are not touched
 elist_c<-backgroundCorrect(elist_c, method="normexp",
 normexp.method="saddle")
 
 #now change it back
 elist2$C<-elist_c$E
 elist2$Cb<-elist_c$Eb
 ##=====>end of background correction
 
 indices<-grep("HumanIgG", elist2$cgenes$Name)
 elist2$cgenes[indices,]
 
#now we need to start checking the error models
# we need to check the block effects

 setwd("E:\\feng\\LAB\\MSI\\AbSynthesis\\proteinArray\\Run2015_10_08\\AnalysisResults");
groups<-c("HumanIgG",
			"HumanIgG", "Anti-HumanIgG", 
			"Anti-HumanIgG1","Anti-HumanIgG2","Anti-HumanIgG3","Anti-HumanIgG4",
			"HumanIgG1","HumanIgG2","HumanIgG3","HumanIgG4"
			)
 #first plot the all IgGs
indices<-grep(groups[1], elist2$cgenes$Name)
expByBlock<-cbind(elist2$C[indices,],elist2$cgenes[indices,c("Block","Name")])

pdf("BlockEffects_100082015_poolAllIgGs.pdf")
op<-par(mfrow=c(3,2),
	pty="m" #"s" for square, "m" for maximal plot area
	)
numArrays<-length(colnames(elist2$C))	
#by each array now
for(i in c(1:length(colnames(elist2$C))))#<-1 ##array one
{
	arrayName<-colnames(elist2$C)[i]
	#colnames(expByBlock)<-paste("a",colnames(expByBlock),sep="")
	boxplot(expByBlock[,i]~Block, data=expByBlock,log="y",
		main=paste("Array ",arrayName, sep=""),
		xlab="Block",ylab="Expression (log)",
		ylim=c(min(expByBlock[,c(1:numArrays)]),max(expByBlock[,c(1:numArrays)])),
		col=i+1
		)
 } 
 par(op)
 dev.off()

 
#second plot the all anti-IgGs
for(k in c(2:length(groups)))
{
	indices<-grep(paste("^",groups[k],sep=""), elist2$cgenes$Name)
	expByBlock<-cbind(elist2$C[indices,],elist2$cgenes[indices,c("Block","Name")])
	numArrays<-length(colnames(elist2$C))	
	
	pdf(paste("BlockEffects_100082015_",groups[k],".pdf",sep=""))
	op<-par(mfrow=c(3,2),
		pty="m" #"s" for square, "m" for maximal plot area
		)
		
	#by each array now
	for(i in c(1:length(colnames(elist2$C))))#<-1 ##array one
	{
		arrayName<-colnames(elist2$C)[i]
		#colnames(expByBlock)<-paste("a",colnames(expByBlock),sep="")
		boxplot(expByBlock[,i]~Block, data=expByBlock,log="y",
			main=paste("Array ",arrayName, sep=""),
			xlab="Block",ylab="Expression (log)",
			ylim=c(min(expByBlock[,c(1:numArrays)]),max(expByBlock[,c(1:numArrays)])),
			col=i+1
			)
	 } 
	 par(op)
	 dev.off()
} 
 
 
 ####code to "debugging"
 ctrData_IgG<-elist3_c$E[indices,]
 indices<-grep("HumanIgG1", elist3_c$genes$Name)
 ctrData_IgG1<-elist3_c$E[indices,]
boxplot(ctrData_IgG1)
 indices<-grep("HumanIgG2", elist3_c$genes$Name)
 ctrData_IgG2<-elist3_c$E[indices,]
boxplot(ctrData_IgG2)
indices<-grep("HumanIgG3", elist3_c$genes$Name)
 ctrData_IgG3<-elist3_c$E[indices,]
 ctrData_IgG3N<-cbind(ctrData_IgG3,elist3_c$genes$Name[indices])
boxplot(ctrData_IgG3)
indices<-grep("HumanIgG4", elist3_c$genes$Name)
 ctrData_IgG4<-elist3_c$E[indices,]
boxplot(ctrData_IgG4)
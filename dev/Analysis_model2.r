
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
elist2$array.type<-"ProtoArray"
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

###calling the PAA to do the normalization
library(PAA) 

elist_norm<-normalizeArrays(elist=elist2, method="rlm",
  cyclicloess.method="fast", control ="HumanIgG")

 #now we have everything, we just need to rewrite the data into a format
 #that can be used by linear regression
 sampleNames<-colnames(elist_norm$E)
 geneNames<-elist_norm$genes[,"Name"]
 
 #now we need to split the data into two pieces, 
 #first 
 dtFrame<-c()
 for(i in 1:length(sampleNames))
 {
	Gene<-geneNames
	Sample<-rep(sampleNames[i],length(geneNames))
	MFI<-elist_norm$E[,i]
	
	temp<-data.frame(Gene,Sample,MFI)
	if(i==1)
	{			
		dtFrame<-temp
	} else
	{
		dtFrame<-rbind(dtFrame, temp)
	}
 }
 
 mdl<-lm(MFI~Gene*Sample, dtFrame)
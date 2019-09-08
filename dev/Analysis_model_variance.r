#testing and showing variations on array

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

 #now we have everything, we just need to summarize the data
 #into "residue" over the sliding window over protein amount
 
 MFI<-elist_norm$E[,1]  #for control now
 proteinAmt<-elist_norm$genes[,"Protein.Amount"]
 
 #now we need to summarize the data to get "residue"
 #in this case, we don't run regression fitting yet,
 #so we only define a sliding window aross protein amount
 #and then average the windowed MFI and get the "residua"
 #first 
 winSize<-50
 maxAmt<-max(proteinAmt)
 numBins<-ceiling(maxAmt/winSize)
 MFI_res<-MFI
 for(i in 1:numBins)
 {
	aveMFI<-ave(MFI[proteinAmt>=(i-1)*winSize&proteinAmt<i*winSize])
	MFI_res[proteinAmt>=(i-1)*winSize&proteinAmt<i*winSize]<-MFI[proteinAmt>=(i-1)*winSize&proteinAmt<i*winSize]-aveMFI
	
 }
 
 #mdl<-lm(MFI~Gene*Sample)#, dtFrame)
 plot(proteinAmt, MFI_res)
 reslm<-lm(MFI_res~proteinAmt)
 abline(reslm,col=2, lwd=2)
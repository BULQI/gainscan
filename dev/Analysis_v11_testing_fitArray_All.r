#code to do analysis of first run of protoArray with isotype control
#Update
#------------
#8/21/2017
#this is the "final" code to run the PMT adjust for processing isotype control data
#1)using the data have been carefully "acquired" with better repeatability. 
#		"AD_Acquire5th.RData" and "AD_BC_Acquire5th.RData" (background corrected data)
#2)with wrapper functions to process them all.
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
library(ARPPA)
pmt_saturated_reading<-65535

#first you need to create a target file including all the meta data for each array, like below
 
 #now try to read the gpr files
 gpr <- "E:\\feng\\LAB\\MSI\\AbSynthesis\\proteinArray\\protoArray\\11_28_2016\\iggisotyprctrl_feng_12052016"
 targets <- list.files(path="E:\\feng\\LAB\\MSI\\AbSynthesis\\proteinArray\\protoArray\\11_28_2016\\iggisotyprctrl_feng_12052016",
					pattern = "^target_second.txt", full.names = TRUE)
#ad.elist <- loadGPR(gpr.path=gpr, targets.path=targets, array.type="ProtoArray",aggregation="none")

#save data
###----> save(ad.elist, file=paste(gpr, "/AD.RData",sep=""), compress="xz")
#load data
 load(paste(gpr,"/AD_Acquire5th.RData",sep=""))
 
 ######==========start doing the 5PL fitting 
 ###visualize the data first to see whether it is a good model
 pmts<-seq(250, 650, by=50)
 #array1.index<-c(1:9)
 
 #get information out
#array1.targets<-ad.elist$targets[array1.index,]

#array1.E<-ad.elist$E[,array1.index]

###let's do background correction first
#ad.elistBC <- backgroundCorrect(ad.elist, method="normexp",
# normexp.method="saddle")
 
 #now we need to rearrange the control and testing proteins to background correct control
 #it works this way that background correction only working on elist$E
 #so we have to do the rearrangement
 #elist_c<-ad.elist
 #elist_c$E<-ad.elist$C
 #elist_c$Eb<-ad.elist$Cb
 #this is good enough, we don't have to change elistC$C, since other fields are not touched
#elist_c<-backgroundCorrect(elist_c, method="normexp",
#normexp.method="saddle")
 
 #now change it back
 #ad.elistBC$C<-elist_c$E
 #ad.elistBC$Cb<-elist_c$Eb
 ##=====>end of background correction
 
#<---section 1	
#----need to fit c, b,g; for fixing both and high ends.
#save data to disk.
 setwd( "E:\\feng\\LAB\\MSI\\AbSynthesis\\proteinArray\\protoArray\\11_28_2016\\iggisotyprctrl_feng_12052016")
#save(ad.elistBC,array1.C,array1.E, file= "AD_BC.RData", compress="xz")
#load data
 load("AD_BC_Acquire5th.RData")
x<-pmts 
#ylog<-F; 

aggregated<-F

#get the testing data by block for testing 
#now getting data ready to do the fitting......
		object<-ad.elistBC
#--		blkNum<-max(object$gene[,"Block"])
#--		arrys<-unique(object$target[,"Array"])
#--		arryNm<-length(arrys)
		
#--		i<-1;j<-1
#--		cat("block ", i,"/",blkNum, "...\n")
		#now take out the array data from the object
#--		aNm<-arrys[j]
		
#--		idx.aNm<-which(object$target[,"Array"]==aNm)
#--		#now we get the index of the one array j, now need to do fitting by block
		#i<-1
		#for the array, we need to get the protein data by block, 4 blocks
		#--idx.prE<-which(is.element(object$gene[,"Block"],c((i):(i+12)))) 
		#idx.prC<-which(is.element(object$cgene[,"Block"],c((i):(i+24))))  
		
		#now get idx, need to get the data into matrix
#--		y<-list()
#--		y$C<-object$C[,idx.aNm] ; #rbind(object$E[idx.prE, idx.aNm], object$C[idx.prC,idx.aNm])
#--		y$E<-object$E[,idx.aNm];
#--		y$cgene<-object$cgene;
#--		y$gene<-object$gene;
		
#--		block.size<-12
		#need to get data in below based on input,
#--		blkNum<-max(y$gene[,"Block"])
#--		if(block.size<1||block.size>blkNum)
#--		{
#--			stop(paste("only a block size between 1 and ", blkNum,sep=""))
#--		}
		
		#prepare the output
#--		fr<-list();#this is the fitting results
#--		dr<-list("E"=y$E[,1], "C"=y$C[,1], "E.adj"=y$E[,1], "C.adj"=y$C[,1], 
#--				"E.gain"=y$E[,1], "C.gain"=y$C[,1]); #this is the gain adjust data output, E/C
															#is the array holding the optimal data for this 
															#array, and PMT.gain is the gain setting
															#for the data
						
#--		j<-1
#--		cat("Fitting data by blocks: ",j, "/",blkNum/block.size, "...\n")
		#for the array, we need to get the protein data by block, 
#--		bg<-(j-1)*block.size+1
#--		ed<-(j-1)*block.size+block.size
#--		if(ed>blkNum)
#--		{
#--			ed<-blkNum
#--		}
		
#--		idx.prC<-which(is.element(y$cgene[,"Block"],c(bg:ed)))
#--		ydata<-y$C[idx.prC,]
#--		fit.mode<-"control"
#--		if(fit.mode=="all")
#--		{
#--			idx.prE<-which(is.element(y$gene[,"Block"],c(bg:ed)))
#--			ydata<-rbind(y$E[idx.prE,], ydata)
#--		}
#--		threshold<-60
#--		data.aggregated<-F
		#now data is ready, rm the low expressed protein, since they are confusing the fitting
#--		ydata5PL<-rmLows(ydata, threshold=threshold, index=which.max(x), aggregated=data.aggregated)
#--		xlen<-length(x)
		
		#now we have the data, reay, aggregated
#--		input<-gainAdjust.prepareInput(ydata5PL,x)
#--		yinput<-input$y
#--		xinput<-input$x
			
		a=1.1;d<-65535
		
		#get the initial values for fitting, by picking the data series within the middle range 
#--		initsLM<-gainAdjust.prepareInitsLM(ydata5PL,x, aggregated=data.aggregated)
		#parStart<-c(6,0.1,0.13,-0.2,-0.5,-0.55, -0.5,1)
		
		#a<-a
		#d<-d
			
		#for three parameters, xmid, scal and g
#--		parStart<-c(6.5,0.005,0.015,initsLM$inits[-initsLM$ref.index]) #<---3pl, keeping a and d constant
#--		data.ylog<-F
		
#--		cat("\tinitial fitting for shifts and parameters.....\n")
#--		gainAdjust.fit<-nls.lm(par=parStart, 
			#fn=gainAdjust.fnRes4pShift,  #<----4pl, varying a
#--			fn=gainAdjust.fnRes3pShift, #<----3pL, keeping a and d constant
			#y=log(yinput), 
#--			y=yinput,
#--			x=log(xinput), xlen=xlen, aggregated=data.aggregated, ylog=data.ylog,a=a,
#--			shift.index=initsLM$ref.index,
#--			model.weight="uniform", order=1, control = nls.lm.control(nprint=1,maxiter=50)
#--		)
#--		cat("\trefining fitting for shifts .....\n")	
		#get the initial values for fitting, by picking the data series within the middle range 
#--		initsLM<-gainAdjust.prepareInitsLM(ydata,x, aggregated=data.aggregated)
		#parStart<-c(6,0.1,0.13,-0.2,-0.5,-0.55, -0.5,1)
#--		input<-gainAdjust.prepareInput(ydata,x, var.log=T)	
#--		fitShift<-fitShifts5plSingle(ydata,x, c(a,d,gainAdjust.fit$par[1:3]), data.aggregated
			#, mvar=input$var
#--			);
		#}
		
#--		aggregated.mode="geometric"
#--		ydata_ag<-gainAdjust.aggregate(ydata, mode=aggregated.mode)
#--		aggregated_3rd<-T
#--		x.aligned<-gainAdjust.alignData(fitShift,x, aggregated=T, model="fix.all", ref.index=initsLM$ref.index)##<---for 4Pl, fixing highest end only

#--		input.aligned<-gainAdjust.prepareInput(ydata_ag,x.aligned)
#--		yinput.aligned<-input.aligned$y
#--		xinput.aligned<-input.aligned$x		
#--		parStartFpl<-c(a,
#--						gainAdjust.fit$par[1],gainAdjust.fit$par[2],gainAdjust.fit$par[3])
#--		cat("\tfinal fitting for parameters.....\n")
#--		gainAdjust.fitFinalP<-nls.lm(par=parStartFpl, 
#--				fn=gainAdjust.fnRes4pFpl,
				#fn=gainAdjust.fnRes3pFpl,
				#y=log(yinput.aligned), 
#--				y=yinput.aligned, 
#--				x=log(xinput.aligned),#xlen=xlen, aggregated=aggregated_3rd,
#--				ylog=ylog, model.weight="sqrt",order=2,d=d,#a=0.1,
#--				control = nls.lm.control(nprint=1, maxit=100)
#--		)	

		
#--		xsim<-seq(4,7,by=0.01)
		
#--		fitShift$par<-gainAdjust.insertAt(fitShift$par,0,initsLM$ref.index)
#--		fitShift$par<-rep(fitShift$par,rep(2, length(fitShift$par)))
		
#--		plot(exp(xsim), f5pl(c(gainAdjust.fitFinalP$par[1],d, gainAdjust.fitFinalP$par[2:4]),xsim), type="l", main="5-p logistic fitting of intensity vs. PMT gain",
#--				xlab="PMT gain (log)",ylab="intensity (log)",log="xy",
#--				lty=1,col=1,lwd=2
#--				)
			#plot(log(xinput.aligned), (yinput.aligned), type="p")
		
#--		lens<-10
#--		for(i in 11:30)
#--		{		
		#lines(x, ydata5PL[i,],col=2)	
		#lines(x, ydata5PL[i,],col=3)
		
#--		lines(exp(log(x)+fitShift$par[i]),ydata[i,],col=i,lty=2)
		#lines(exp(log(x)+fitShift$par[958/2-1]),ydata5PL[906,],col=3,lty=2)
#--		}
		
		
		
		
		
		
		
#--		system.time(rs<-gainAdjust.fit5plArray(y, x, data.ylog=F, data.aggregated=F,
#--		fit.mode="control",
#--		fit.aggregated=F, aggregated.mode="geometric", order=1,
#--		block.size=12,a=1.1, d=65535,threshold=60,
#--		debug=T
#--		))
		
		###now call function to do the fitting for 
		system.time(
		yAdj<-PMT.adjust(object, #raw date after reading the GPR files
			x, #vector holding the PMT gains for reading the chips, the arrangement of data in PMT reading is identical to the element order in x.
			data.ylog=F, data.aggregated=F,#description of data, usually, don't need to change, for advance user only
			fit.mode="control", ##fit the 5pl with control data only, to speed up 
			fit.aggregated=F, #fit the 5pl with aggregated data, for advanced data
			aggregated.mode="geometric", #controlling the way to do aggregation, geometric is preferred
			block.size=12, #size of blocks to do fitting, the more the better??. but too slow, 12 is appropriate
			a=1.1, d=65535, #change only when you know what you do.
			debug=T  #showing the debugging information
		)
		);







		
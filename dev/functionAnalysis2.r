#function to process the simulated data, 
# feng@BU 08/30/2017
#
# 0) revise the model Y=a+b+g+k*e    scale and location model, effects of both additive and multiplicative
#
#1) read the data from the disk
#1a) estimate the overall variance as a surgate variance of base level
#1b) estimate k by Empirical bayes with a sample variance versus population variance (following a ki-square distribution)
#1c) derive posterior distribution of k
#1d) rescale/transform y by scaled off k to lead to equal variance model 
#2) process to see the linear regress results
#		pick the top portions to see the identification performance
#		fdr
#		roc
#		
#3) do the ordinary t test to compare the performance of linear regression model
#4) sumarize the results
#5) plotting and reporting
###############################
library(ARPPA)
library(qvalue)
#install.packages("pROC")
library(pROC)


#filename is the base name of the file. the file should be 
#		r data file (".RData"). filename="result_"
#sample.size, used for t test, the sample size for each group
#percent.nonZero, indicate the proportion of nonzero interaction in the simulation 
#object.load, used to indicate the loaded object name
#mode, 1 for comparison with either negative control or isotype control, other values for no control
analyzeData<-function(path,filename="result_", repeats=1, 
					sample.size=5,
					proportion.nonZero=0.02,
					object.load="lstToSave",
					mode=1
			)
			
{
#path<-filePath
sampleSize=sample.size;
prob.nonZero=proportion.nonZero;
cat("start doing the data reading......\n");
flush.console();
repeats<-repeats
setwd(path)
#("E:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\EqualVar_NegativeCon");
#setwd("~/Desktop/arppr/EqualVar_NegativeCon/")
i<-19
portions<-seq(0.000, 1.0, 0.001)
stats_df<-data.frame("portions"=portions);
stats_qval_tpr<-data.frame("cutoff"=portions)
stats_qval_fpr<-data.frame("cutoff"=portions)

stats_qval_tpr_ordT<-data.frame("cutoff"=portions)
stats_qval_fpr_ordT<-data.frame("cutoff"=portions)

roc.auc<-rep(0,repeats)
roc_ordT.auc<-rep(0,repeats)
for(i in c(1:repeats))
{
	
	cat("reading the ",i, "/", repeats," data sets.....\n")
	#load(paste(filename, i,"_rlm.RData",sep=""))
	load(paste(filename, i,".RData",sep=""))
	
	cat("count the correct ones.....\n");
	flush.console();
	#depending on which data,
	lstData<-get(object.load);#lstToSave
	eval(paste("rm(",object.load,")",sep=""))
	
	#estimate k posterior------------>
	dataExp_EN<-dataExp_EN_List[[1]]
	sMean<-sampleMean(dataExp_EN, nTreatment, sampleSize)
	sVar<-sampleVariance(dataExp_EN, nTreatment, sampleSize)
	sdf<-sampleSizes(dataExp_EN, nTreatment, sampleSize)-1
	#get the sample variance prior by assuming a scale-chisquare
	s2<-as.vector(as.matrix(sVar))
	df<-as.vector(as.matrix(sdf))
	
	r2<-sampleVariancePrior(s2,df)
	s02<-r2[[2]]
	df0<-r2[[1]]
	
	#calculate the posterior variance
	var_pos<-sampleVariancePosterior(sdf,sVar,  df0, s02 )
	var_pos_vec<-(as.vector(t(var_pos)))
	var_pos_vec_all<-rep(var_pos_vec,rep(5, length(var_pos_vec)))
	
	-->>>>>>>
	#try to know which one is nonZero as the true values
	gammaTrue<-lstData$data$gamma
	if(mode==1)
	{
		gammaTrueIndex<-which(gammaTrue[,2]!=0)
	} else
	{
		gammaTrueIndex<-which(gammaTrue[,2]!=0|gammaTrue[,1]!=0)
	}
	
	#now need to sort the p-value of the 
	#and pick top n genes to be the significant ones as the analysis results
	#
	lr_coe<-lstData$lregCoef
	lr_coe_names<-rownames(lr_coe)
	lr_coe_gamma_index<-which(grepl("gene\\d*:group2",lr_coe_names)); #indices for interaction terms only
	lr_coe_int<-lr_coe[lr_coe_gamma_index,];#get entries for interaction terms
	lr_coe_int_names<-lr_coe_names[lr_coe_gamma_index];

	#sort the array according to the prob
	lr_coe_int_prob<-lr_coe_int[,"Pr(>|t|)"]
	#lr_coe_int_prob<-(1-pt(abs(lr_coe_int[,"t value"]),df=2000*2*5-4000))*2
	
	lr_coe_int_sort_order<-order(lr_coe_int_prob) #NOTE:the names is one more than its index, eg. 811 is gene812:group2
	#lr_coe_int_sort<-lr_coe_int[lr_coe_int_sort_order,]
	
	#lr_coe_int_prob_qvalue<-qvalue(lr_coe_int_prob[lr_coe_int_sort_order])$qvalue
	lr_coe_int_prob_qvalue<-tryCatch(
		qvalue(lr_coe_int_prob)$qvalue,   #(lr_coe_int_prob[lr_coe_int_sort_order])$qvalue,
		error=function(c){
		cat("******ERROR in calling \"qvalue\"\n")
		cat("\t",c$message,"\n")
		cat("calling p.adjust instead\n")
		lr_coe_int_prob_qvalue<- p.adjust(lr_coe_int_prob, method="BH")#(lr_coe_int_prob[lr_coe_int_sort_order], method="BH")		
		},
		warning = function(c) c$message,
		message = function(c) c$message
	)
	
	roc_response<-rep(0, length(gammaTrue[,1]))
	roc_response[gammaTrueIndex]<-1
	
	roc_predict<-lr_coe_int_prob_qvalue #qvalue(lr_coe_int_prob)$qvalue
	roc.obj<-roc(roc_response[-1], roc_predict)
	roc.auc[i]<-roc.obj$auc
	

	#now we have everything, just need to collect statistics
	p_correct<-portions;#initialize the vector
	p_correct_byQ<-portions;#initialize the vector
	p_fpr_byQ<-portions;#initilize
	
	p_correct_byQ_ordT<-portions;#initialize the vector
	p_fpr_byQ_ordT<-portions;#initialize
	
	
	#now we are doing -->ordinary t tests<--
	#get the data first
	dtExp<-lstData$data$exp
	prob_ordT<-rep(0,length(dtExp[,1]))
	#run t.test
	for(nn in c(1:length(dtExp[,1])))
	{
		prob_ordT[nn]<-t.test(dtExp[nn,c(1:sampleSize)],dtExp[nn,c((1+sampleSize):(sampleSize+sampleSize))])$p.value
	}
	prob_ordT_Q<-	tryCatch(
		qvalue(prob_ordT)$qvalue,
		error=function(c){
		cat("******ERROR in calling \"qvalue\"\n")
		cat("\t",c$message,"\n")
		cat("calling p.adjust instead\n")
		prob_ordT_Q<- p.adjust(prob_ordT, method="BH")		
		},
		warning = function(c) c$message,
		message = function(c) c$message
	)
	roc_ordT<-roc(roc_response, prob_ordT_Q)
	
	roc_ordT.auc[i]<-roc_ordT$auc
	j<-1
	for(j in c(1:length(portions)))
	{
		if(mode==1)
			factors<-1
		else
			factors<-2
		#for each portion, we need to check what is percentage to be correct
		numToPick<-floor(lstData$data$params[1]*portions[j]*factors)
		#picking from the order array
		genesToPick_byProb<-lr_coe_int_sort_order[c(1:numToPick)]+1 #see NOTE above, index is one more less than the name
		#now we got the genes, just need to check whether the genes so far are the TRUE ones
		numOfCorrect<-sum(is.element(genesToPick_byProb, gammaTrueIndex));
		p_correct[j]<-numOfCorrect/numToPick
		
		#now doing the q value, using the portion as the qvalue cutoff
		numToPick_byQ<-floor(portions[j]*lstData$data$params[1]*factors);#sum(lr_coe_int_prob_qvalue<=portions[j])
		lr_coe_int_Q_sort_order<-order(lr_coe_int_prob_qvalue)	
		genesToPick_byQ<-lr_coe_int_Q_sort_order[c(1:numToPick_byQ)]+1
		#now we got the genes, just need to check whether the genes so far are the TRUE ones
		numOfCorrect_byQ<-sum(is.element(genesToPick_byQ, gammaTrueIndex));
		p_correct_byQ[j]<-numOfCorrect_byQ/floor(lstData$data$params[1]*lstData$data$params[3]*factors)
		p_fpr_byQ[j]<-(length(genesToPick_byQ)-numOfCorrect_byQ)/floor(lstData$data$params[1]*(1-lstData$data$params[3]*factors))
		
		#now doing ordinary T tests<--------
		numToPick_byQ_ordT<-numToPick_byQ;#sum(prob_ordT_Q<=portions[j])
		prob_ordT_Q_sort<-order(prob_ordT_Q)
		geneToPick_byQ_OrdT<-prob_ordT_Q_sort[c(1:numToPick_byQ_ordT)]#this is different, there is +1 here!!!!!!!
		numOfCorrect_byQ_ordT<-sum(is.element(geneToPick_byQ_OrdT, gammaTrueIndex));
		
		p_correct_byQ_ordT[j]<-numOfCorrect_byQ_ordT/floor(lstData$data$params[1]*lstData$data$params[3]*factors)
		p_fpr_byQ_ordT[j]<-(length(geneToPick_byQ_OrdT)-numOfCorrect_byQ_ordT)/floor(lstData$data$params[1]*(1-lstData$data$params[3]*factors))
	}
	stats_df[,paste("repeat_",i, sep="")]<-p_correct;
	stats_qval_tpr[,paste("repeat_",i,sep="")]<-p_correct_byQ;
	stats_qval_fpr[,paste("repeat_",i,sep="")]<-p_fpr_byQ;
	
	stats_qval_tpr_ordT[,paste("repeat_",i,sep="")]<-p_correct_byQ_ordT;
	stats_qval_fpr_ordT[,paste("repeat_",i,sep="")]<-p_fpr_byQ_ordT;
	cat("done!\n")
}
cat("start doing the plotting and reporting...\n")
flush.console();
#statistics with the stats_df array
mean_correct<-portions
max_correct<-portions
min_correct<-portions
std_correct<-portions
m_stats_df<-as.matrix(stats_df)
for(k in c(1:length(portions)))
{
	mean_correct[k]<-mean(m_stats_df[k,c(-1)])
	max_correct[k]<-max(m_stats_df[k,c(-1)])
	min_correct[k]<-min(m_stats_df[k,c(-1)])
	std_correct[k]<-sqrt(var(m_stats_df[k,c(-1)]))
}

stats_df[,"mean"]<-mean_correct
stats_df[,"min"]<-min_correct
stats_df[,"max"]<-max_correct
stats_df[,"std"]<-std_correct
#op<-par(mfrow=c(2,1),
#	pty="s" #"s" for square, "m" for maximal plot area
#	)
plot(c(0.001,0.51),c(0.02, 1.1), type="n", main="true discovery rate for protein selection", 
	xlab="portion of protein selected", ylab="portion of true discovery", log="x")
lines(stats_df[,1], stats_df[,"mean"], col=1, lty=1, lwd=2)
lines(stats_df[,1], stats_df[,"min"], col="grey", lty=2, lwd=1)
lines(stats_df[,1], stats_df[,"max"], col="grey", lty=2, lwd=1)
lines(c(prob.nonZero, prob.nonZero), c(0,1.2), col=2, lty=3)

#show stats
stats_df[c(1:30),c("portions", "mean", "min", "max", "std")]
############adding code to do FDR

#statistics with the roc array
mean_tpr<-portions
mean_fpr<-portions
for(k in c(1:length(portions)))
{
	mean_tpr[k]<-median(as.matrix(stats_qval_tpr[k,c(-1)]))
	mean_fpr[k]<-median(as.matrix(stats_qval_fpr[k,c(-1)]))
}


plot(c(0,1), c(0,1),type="n", main="C: No Control", ylab="True Positive Rate", xlab="False Positive Rate")
#for(ii in c(1:repeats))
#{
#	lines( stats_qval_fpr[,ii+1],stats_qval_tpr[,ii+1],col="red", lty="dotted")
#}
lines(mean_fpr,mean_tpr, col="red",lty=1,lwd=2)
#legend(0.7,0.2,c("mean ROC"),col=c("red"), lty=c(2),lwd=c(2))
#text(0.7,0.3,labels=paste("AUC:",mean(roc.auc)))

#statistics with the roc array for ---->ordT
mean_tpr_ordT<-portions
mean_fpr_ordT<-portions
for(k in c(1:length(portions)))
{
	mean_tpr_ordT[k]<-median(as.matrix(stats_qval_tpr_ordT[k,c(-1)]))
	mean_fpr_ordT[k]<-median(as.matrix(stats_qval_fpr_ordT[k,c(-1)]))
}


#plot(c(0,1), c(0,1),type="n", main="ROC", ylab="True Positive Rate", xlab="False Positive Rate")
#for(ii in c(1:repeats))
#{
#	lines( stats_qval_fpr_ordT[,ii+1],stats_qval_tpr_ordT[,ii+1],col="light green", lty="dotted")
#}
lines(mean_fpr_ordT,mean_tpr_ordT, col="green",lty=1,lwd=2)
legend(0.7,0.2,c("model", "ordinary t"),col=c("red", "green"), lty=c(2, 2),lwd=c(2,2))
text(0.7,0.25,labels=paste("AUC ordinary t:",mean(roc_ordT.auc)))
text(0.7,0.35,labels=paste("AUC model:",mean(roc.auc)))
#par(op);

####starting here, prepare the output
	
	retDataList<-list("roc.auc.model"=roc.auc,"roc.auc.t"=roc_ordT.auc, 
				"qvalule.tpr.model"=stats_qval_tpr, "qvalule.fpr.model"=stats_qval_fpr,
				"qvalule.tpr.t"=stats_qval_tpr_ordT, "qvalule.fpr.t"=stats_qval_fpr_ordT,
				"fdr"=stats_df 
				#, "param.alpha"=lstData$data$alpha,"param.beta"=lstData$data$beta,
				#"param.gamma"=lstData$data$gamma,
				);
	return(retDataList)
}#end of function

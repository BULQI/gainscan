###this is the module to define the data simulation and related functions for expression
#####this is r developing code for running simulation to generate 
##the DNA/protein microarray data
## the model (ref: smyth 2004, Bayesian hierarchical model and linear model)
##  linear model:
##             Y_ijk=alpha_i+beta_j+ gamma_ij+ epsilon_ijk
##		see the doc named: "Microarray data simulation.doc"
##    feng@BU   09/03/2016
####################################

#note: all the effects are fixed, although the observed ones are "random", since there is 
# 	 observation error or disturbance. the error or disturbance made the estimated effects 
#    "random". 
##S3 function

#comment for Oxygen2 helper page
#'@title S3 function to simulate expression data
#'@description \code{simulateExpression} simulates the expression data assuming
#'	a specific/non-specific effect model (biologically) and a Bayesian 
#'  hierarchical model (statistically)
#'@details This function is used to simulate the data assuming statistically a 
#'  a linear model with Bayesian and hierarchical structure according to Smyth 2004.
#'  Biologically, it parses the observed expression data into effects from different 
#'	sources, mainly nonspecific and specific ones. The nonspecific effects comes 
#'  due to the nonspecifi binding by either proteins or treatment agents. The specific one 
#'  is cause by the specific interaction due to protein and treatment agents (here antibodies).
#'  Linear model:  >>>>>>>>>To be continue>>>>>>>>>>>>>>
#'	About the control groups in the simulation, it is always true that there are 
#'	gene-wise control groups. Under the real condition, these are the genes that are
#'	printed at very little amount on the array. In case of the controls in the treatment group,
#'	there are 3 cases. First is the isotype control group, in which the isotype 
#'	antibody shows no reactivity to genes, but only nonspecific effects. Second 
#'	one is the negative control, meaning there is no antibody at all but only solution.
#'	in this case, there is no nonspecific or interaction effects. The third case
#'	has no control. "Therefore, we need to compare the every treatment group with
#'	the grand mean to estimate the interaction. The assumption here is that 
#'	non-zero interaction is a rare event and grand mean should be a good estimation
#'	of the nonspecific effects of the gene."
#'@param nGene numeric the number of genes to be simulated
#'@param nTreatment number of treatment groups, eg number of antibodies
#'@param sampleSize number of replicates for each gene by treatment group
#'@param control.isotype logic indicating whether for the treatment group
#'	there is isotype controls.
#'@param control.negative logic indicating whether for the treatment group
#'	there is negative controls, meaning the vehicle control and no antibody.
#'	control.negative and control.isotype are exclusive. if both control.negative
#'	and control.isotype are set to be true, then the negative case is assumed.
#'	When alpha and beta are set, all the control arguments are ignored.
#'	Again, there are always gene control groups and we always set first one
#'	to be the control group by gene.
#'
#'@param alpha numeric vector specifying the non-specific gene effects with 
#'  length of nGene
#'@param alpha.mean and alpha.sigma numeric the parameters used to randomly specify
#	the non-specific alpha effects following the normal distribution with
#	mean and standard deviation. When parameter alpha is specified by
#	the user. These two parameters are ignored
#'@param beta numeric vector specifying the non-specific treatment/antibody
#'  effects with length of nTreatment
#'@param beta.mean and beta.sigma numeric the parameters used to randomly specify
#	the non-specific beta effects following the normal distribution with
#	mean and standard deviation. When parameter alpha is specified by
#	the user. These two parameters are ignored
#'@param gamma numeric The specific effect between gene and treatment.
#'	Under the null hypothesis, this is zero assuming no interacting effects.
#' 	The observed expressions for many of the protein and treatment are
#'	caused mainly because of non-specific effects with disturbance. 
#'@param prob.nonZero numeric this is true partion for those proteins
#' 	having the non-specific interaction with the treatment agents.
#'@param gamma.sigma numeric the parameter specific the distribution from which
#'	the non-zero gamma drawing values. Here we assuming a Gaussian 
#'	distribution with mean of zero and standard deviation of 
#'	gamm.sigma. If gamma is specified by the user, then gamma.sigma and 
#'	prob.nonZero will be ignored.
#'@param epsilon.si numeric matrix the standard deviations for Gaussian
#'	distributed errors/disturbances.	
#'@param epsilon.d0 and epsilon.s0 numeric The hyperparameter specifying
#'	the prior distribution of for the variance for each inidividual group.
#'	Here we assume a heteroscedastical error model. 
#'		episilon ~ N(0, si2)
#'		d0*s02/si2~chiSquare(d0,di)
#'		Yijk=Yijk_hat+epsilon 
#'		di is the df of si2
#'
#'@return a list containing 1)a numeric data matrix with a format of gene by treatment/antibody. 
#'	The columns are arrangend to have repeats; 2)parameters alpha for genes;
#'	3)parameters beta for treatment; 4)parameters gamma for interaction;
#'	5)vector of other parameters (named), number of Genes, number of Treatment/Group, 
#'	probability for nonzero interactions;
#'	6)vector of data variance; 7)sample.size (balanced or unbalanced design)
#'	 
#'@examples 
#' 
#' nGene<-100
#' nTreatment<-2   #number of different beta
#' sampleSize<-5 ##this is the number of repeats for each group

#' alpha.mean<-0  #variance for alpha prior
#' alpha.sigma<-3
#' beta.mean<-0
#' beta.sigma<-2   #variance for beta prior
#'
#' gammaN0.sigma<-10  # the unscaled factor for gamma given gamma<=>0
#' p_diff<-0.01

#' #priors for variance distribution
#' d0<-5
#' s0<-2

#' #call it
#' set.seed(2004);
#' dataExp<-simulateExpression(nGene, nTreatment, sampleSize,
#'					alpha.mean=alpha.mean, alpha.sigma=alpha.sigma,
#'					beta.mean=beta.mean, beta.sigma=beta.sigma,
#'					prob.nonZero=0.01, gamma.sigma=gammaN0.sigma,
#'					epsilon.d0=d0, epsilon.s0=s0
#'					)
#'
#'@export

simulateExpression<-function(nGene, nTreatment, sampleSize=2,
			alpha=NULL, alpha.mean=NULL, alpha.sigma=NULL, 
			control.isotype=TRUE, control.negative=TRUE, 
			control.index=c(1),
			beta=NULL,beta.mean=NULL, beta.sigma=NULL, 
			gamma=NULL, prob.nonZero=NULL, gamma.sigma=NULL,
			epsilon.si=NULL, epsilon.d0=NULL, epsilon.s0=NULL
			)
			{
				#first need to check for the correct input
				#assuming the nGene and nTreatment and sampleSize
				#set.seed(2004)
				if(sampleSize<2)
				{
					stop("ERROR: sample size has to be bigger than 2!");
				}
				if(nTreatment<1)
				{
					stop("ERROR: number of treatment groups has to be bigger than 1!");
				}
				if(nGene<1)
				{
					stop("ERROR: number of genes has to be bigger than 1!");
				}
				if(nGene<100)
				{
					cat("WARNING: number of gene is small. The estimation might be inconsistent!\n")
				}
				
				if(is.null(alpha))
				{
					if(is.null(alpha.mean)||is.null(alpha.sigma))
					{
						stop("ERROR: no values has been set for alpha (the non-specific effects for genes!");
					}
					#now set alpha according to the distribution
					alpha<-rnorm(nGene,alpha.mean, alpha.sigma)
					alpha[1]<-0;
				}
				if(length(alpha)==1)
				{
					alpha<-rep(alpha, nGene)
				}
				if(length(alpha)!=nGene)
				{
					stop("ERROR:the length of input argument 'alpha' is not correct.")
				}
				#cat("alpha:", alpha,"\n")
				if(is.null(beta))
				{
					if(is.null(beta.mean)||is.null(beta.sigma))
					{
						stop("ERROR: no values has been set for beta (the non-specific effects for treatment group!");
					}
					#now set alpha according to the distribution
					beta<-rnorm(nTreatment,beta.mean, beta.sigma)
					if(control.negative)
					{
					#cat("yes\n")
						beta[control.index]=0;
					}
					if(control.isotype)
					{
					#cat("yes2\n")
						#beta[control.index]
						#here we don't have to do anything.
					}
				}
				if(length(beta)==1)
				{
					beta<-rep(beta,nTreatment);
				}
				if(length(beta)!=nTreatment)
				{
					stop("ERROR:the length of input argument 'beta' is not correct.")
				}
				#cat("beta:",beta,"\n")
				if(is.null(gamma))
				{
					if(is.null(prob.nonZero)||is.null(gamma.sigma))
					{
						stop("ERROR: no values has been set for alpha (the non-specific effects for genes!");
					}
					#now set alpha according to the distribution
					gamma<-matrix(rep(0,nGene*nTreatment),nrow=nGene, byrow=T)
					##in this following code snip, we randomly distribute the differential gamma into
					##different positions across different treatment group
					##other part is 
					numGene_diff<-floor(nGene*prob.nonZero)
					#cat("numGene_diff:",numGene_diff,"\n")
					antibodyIndex<-c(1:nTreatment)
					if(control.isotype||control.negative)
					{
						
						antibodyIndex<-antibodyIndex[-1*control.index]
					}
					cat("antibodyIndex:", antibodyIndex, "\n");
					for(i in antibodyIndex)
					{
						pos_diff<-sample(c(2:nGene), size=numGene_diff, replace=F)
						#cat("pos_diff:",pos_diff,"\n")
						#cat("i:",i,"\n")
						#cat("gamma sub:", gamma[pos_diff,i],"\n")
						#we want to be replace false, since otherwise we will have the collision (same number drawn multiple time)
						gamma[pos_diff,i]<-(rnorm(numGene_diff,0,gamma.sigma))
					}
					#gamma[0,]<-0;
					#if(control.isotype||control.negative)
					#{
					#	gamma[,index]<-0;
					#}
				}
				#cat("gamma:\n")
				#print(gamma)
				#cat("\n")
				if(length(gamma)==1)
				{
					gamma <-matrix(rep(gamma, nGene*nTreatment), nrow=nGene, byrow=T);
				}
				if(dim(gamma)[1]!=nGene||dim(gamma)[2]!=nTreatment)
				{
					stop("ERROR:the length of input argument 'gamma' is not correct.")
				}
				if(is.null(epsilon.si))
				{
					if(is.null(epsilon.d0)||is.null(epsilon.s0))
					{
						stop("ERROR: no distribution parameters has been set!");
					}
					##now ready to generate variance
					epsilon.si<-rScaledChisq(nGene*nTreatment,epsilon.d0,epsilon.s0*epsilon.s0)
					#epsilon.si<-matrix(rchisq(nGene*nTreatment, df=epsilon.d0),nrow=nGene, byrow=T)
					#epsilon.si<-1/(epsilon.d0*epsilon.s0*epsilon.s0)*epsilon.si
					#epsilon.si<-1/epsilon.si
					epsilon.si<-sqrt(epsilon.si)
					epsilon.si<-matrix(epsilon.si,nrow=nGene, byrow=T)
				
				}
				if(length(epsilon.si)==1)
				{
					epsilon.si<-matrix(rep(epsilon.si,nGene*nTreatment), nrow=nGene, byrow=T);
				}
				if(dim(epsilon.si)[1]!=nGene||dim(epsilon.si)[2]!=nTreatment)
				{
					stop("ERROR:the length of input argument 'epsilon.si' is not correct.")
				}
				##with the variance we can generate the data
				#now we have everything ready, do the observation values
				#put them into matrix first
				groupInfo<-rep(0,sampleSize*nTreatment)
				Y_ijk<-matrix(rep(0,nGene*nTreatment*sampleSize),nrow=nGene,byrow=T)
				for(j in c(1:nTreatment)) #for different treatment
				{
					samplePos<-c(1:sampleSize)+(j-1)*sampleSize;
					groupInfo[samplePos]<-j;
					for(k in c(1:nGene))
					{
						#write the meta data first
						Y_ijk[k,samplePos]<-alpha[k]+beta[j]+gamma[k,j]
						Y_ijk[k,samplePos]<-Y_ijk[k,samplePos]+rnorm(sampleSize,0,epsilon.si[k,j])
					}
				}

				#for now, return the matrix and will try other later
				retList<-list( "exp"=Y_ijk, "alpha"=alpha, "beta"=beta, 
						"gamma"=gamma, 
						"params"=c("nGene"=nGene, "nGroup"=nTreatment, "prob.nonZero"=prob.nonZero),
						"sample.Size"=sampleSize,  "epsilon.sd"=epsilon.si
						);
			}
			
#accessary function for converting the matrix to data.frame
#
#Oxygen 2 comment
#'@title convert a data matrix into data frame
#'@description takes in data matrix and convert into a data frame with necessary 
#'	categorical information
#'@details It assumes a data matrix with a format of gene by groups
#'	and then turn it into a data frame with data fields of expression followed
#'	by group id and gene id.
#'@param x matrix containing the data with a format of genes as row and group as column
#'@param nGroup numeric the number of groups
#'@param sampleSize vector the number of repeats in each group. If it is a balanced 
#'	design, then a scalar value is enough
#'@return a data frame holding the data and categorical information
#'@export
matrix2dframe<-function(x, nGroup, sampleSize)
{
	if(missing(x))
	{
		stop("ERROR:undefined function argument \'x\'!")
	}
	if(missing(nGroup))
	{
		stop("ERROR:undefined function argument \'nGroup\'!")
	}
	if(missing(sampleSize))
	{
		stop("ERROR:undefined function argument \'sampleSize\'!")
	}
	
	#now check for the input
	if(!is.matrix(x)&&!is.data.fram(x))
	{
		#not a matrix
		stop("ERROR: the input data is not in a correct fromat!!")
	}
	#make sure the 
	if(length(sampleSize)==1)
	{
		sampleSize<-rep(sampleSize,nGroup)
	}
	#now check the size of group
	if(dim(x)[2]!=sum(sampleSize))
	{
		#something wrong
		stop("ERROR:the input data object doesn't have the correct dimension as specified!.")
	}
	#now , so far so good
	#do the conversion column by column
	dtm<-data.frame()
	gene<-c(1:dim(x)[1])
	index<-0
	for(i in 1:nGroup)
	{
		for(j in 1:sampleSize[i])
		{
			index<-index+1
			exp<-x[,index]
			group<-i
			temp<-data.frame(exp,group, gene)
			if(index==1)
			{
				dtm<-temp
			}
			else
			{
				dtm<-rbind(dtm, temp)
			}
		}
	}
	dtm$gene<-as.factor(dtm$gene)
	dtm$group<-as.factor(dtm$group)
	cat("Done!\n");
	return(dtm);
}

######this is the function used to generate a random sample following gaussian distribution
##### of gamma parameters (the interaction term in the model for protein array)
#

generateGammaMatrix<-function(nGene, nTreatment, prob.nonZero=0,group.differential=NULL)
{
	#now set alpha according to the distribution
	gamma<-matrix(rep(0,nGene*nTreatment),nrow=nGene, byrow=T)
	##in this following code snip, we randomly distribute the differential gamma into
	##different positions across different treatment group
	##other part is 
	numGene_diff<-floor(nGene*prob.nonZero)
	for(i in c(1:nTreatment))
	{
		pos_diff<-sample(nGene, size=numGene_diff, replace=T)
		gamma[pos_diff,i]<-(rnorm(numGene_diff,0,gamma.sigma))
	}
	return(gamma)
}

########function to do some accessory function of average, log transformation
# 
#input
#  object, EList or data matrix. for data matrix, we simply
#  log, if TRUE, take the logarithm of the data first before doing aggregate
#  method, be careful here, if take the logarthm of the data, then no need to run Geometric mean
#	when the input is EList we don't do anything on background table, assuming the background 
#	correction has been done. 
#'@export
array.aggregate<-function(
		object, log=TRUE, log.base=exp(1),
		array.type=c("ProtoArray","unknown"),
		method=c("Geometric","Arithmetic", "None")
	)
	{
		#check the input
		if(is.null(object))
		{
			stop("ERROR in the input, object can not be null")
		}
		if(class(object)!="EListRaw"&&class(object)!="matrix")
		{
			stop("ERROR: unknown object data, can not proceed!")
		}
		if(!is.logical(log))
		{
			stop("ERROR: log must be logical, TRUE or FALSE")
		}
		#log<-match.arg(log)
		array.type<-match.arg(array.type)
		method<-match.arg(method)
		
		if(method=="Geometric"&&log)
		{
			cat("****WARNING: calling to do geometric mean on log-transformed data, NaN could be generated!!\n")
		}
		
		#cat("here before switch\n")
		switch(class(object),
		EListRaw={
				#cat("within ELIST before calling\n")
				object$E<-matrix.aggregate(object$E,log, log.base,method);
				object$genes<-dataframe.aggregateGenes(object$genes);
				if(array.type=="ProtoArray"||object$array.type=="ProtoArray")
				{
					object$C<-matrix.aggregate(object$C,log,log.base, method);
					object$cgenes<-dataframe.aggregateGenes(object$cgenes);
				}
			},
		matrix={#in here, we simply do the aggregation
				#cat("matrix in switch\n")
				object<-matrix.aggregate(object, log,log.base, method)
			}
		)
		
		object
	
	}

#accessory function aiming to do the real aggregation on data matrix
#assuming there are duplicated data next to each other.
#if choose to do with method=None, it simply run log transformation	
	matrix.aggregate<-function(m, log=TRUE, log.base=exp(1),
			method=c("Geometric", "Arithmetic", "None")
		)
	{
		#check the input
		if(is.null(m))
		{
			stop("ERROR in the input, object can not be null")
		}
		if(class(m)!="matrix")
		{
			stop("ERROR: input is not a matrix object, can not proceed!")
		}
		
		#
		if(!is.logical(log))
		{
			stop("ERROR: log must be logical, TRUE or FALSE")
		}
		#log<-match.arg(log)
		#array.type<-match.arg(array.type)
		method<-match.arg(method)
		
		if(method=="Geometric"&&log)
		{
			cat("****WARNING: calling to do geometric mean on log-transformed data, NaN could be generated!!\n")
		}
		
		if(log)
		{
			m<-log(m, base=log.base)
		}
		switch(method,
			Geometric={
					m<-sqrt(m[seq(1,dim(m)[1],by=2),]*m[seq(2,dim(m)[1],by=2),])
			},
			Arithmetic={
				m<-(m[seq(1,dim(m)[1],by=2),]+m[seq(2,dim(m)[1],by=2),])/2
			},
			None={
			}
		)
		m
	}#end of function
	#function to "aggregate" the gene description table of EList object in order to 
	#make the two table referencable
	dataframe.aggregateGenes<-function(dtf)
	{
		#check the input
		if(is.null(dtf))
		{
			stop("ERROR in the input, object can not be null")
		}
		dtf<-dtf[seq(1,dim(dtf)[1],by=2),]
		
	}
	
		#rewrite the data into a R style of data frame by combining the data and info from 
		#matrix and gene data frame
	#input:
	#	datamatrix, datamatrix holding the signals for proteins/genes
	#			rows are genes and columns are arrays
	#			also assuming the matrix has been correctly aggregated
	#			most likely this is the control protein fluorescence matrix
	#	geneframe, data frame holding the genes information
	#			it is also aggregated together with datamatrix
	#	col.indices, this is the indices for arrays to be included in
	#			the analysis. If missing, assuming all columns to be used
	####  '   @return the data frame holding the CVs for the control 
	#'@export
	combineDataNGene<-function(datamatrix, geneframe, col.indices
	)
	{
		if(missing(datamatrix))
		{
			stop("please specify the input data frame.")
		}
		if(missing(geneframe))
		{
			stop("please specify the input gene information frame.")
		}
		if(missing(col.indices))
		{#means doing all
			col.indices<-c(1:length(dataframe[1,]))
		}
		#start doing the variance intra-array with controls
		#intra-array variance first
		sboner_cDtf<-datamatrix
		sboner_cgDtf<-geneframe
		
		#now for intra-array we pick PMT=700 to calcuate the intra-array variance, 
		#since PMT=600 is recommended by the manufacturer.
		#cn<-colnames(sboner_cDtf)
		sboner_pmtHiIndex<-col.indices
		
		sboner_cdat<-data.frame()
		
		#do summarizing and plotting
		for(i in sboner_pmtHiIndex)
		{
			
			#make the raw data ready into a dataframe 
			#so to plot
			sboner_cdat_dtf<-cbind(F=sboner_cDtf[,i],sboner_cgDtf,array=i)
			#edat_dtf<-cbind(F=eDtf[,i],egDtf,array=i)
			
			if(length(sboner_cdat)==0)
			{
				sboner_cdat<-sboner_cdat_dtf
			}else
			{
				sboner_cdat<-rbind(sboner_cdat, sboner_cdat_dtf)
			}
		}#end of for loop
		
		#now return 
		sboner_cdat
	}
	#calculate intra-array cvs
	#input:
	#	datamatrix, datamatrix holding the signals for proteins/genes
	#			rows are genes and columns are arrays
	#			also assuming the matrix has been correctly aggregated
	#			most likely this is the control protein fluorescence matrix
	#	geneframe, data frame holding the genes information
	#			it is also aggregated together with datamatrix
	#	col.indices, this is the indices for arrays to be included in
	#			the analysis. If missing, assuming all columns to be used
	#### '   @return the data frame holding the CVs for the control 
	#'@export
	intraArrayCVs<-function(datamatrix, geneframe, col.indices
	)
	{
		if(missing(datamatrix))
		{
			stop("please specify the input data frame.")
		}
		if(missing(geneframe))
		{
			stop("please specify the input gene information frame.")
		}
		if(missing(col.indices))
		{#means doing all
			col.indices<-c(1:length(dataframe[1,]))
		}
		#start doing the variance intra-array with controls
		#intra-array variance first
		sboner_cDtf<-datamatrix
		sboner_cgDtf<-geneframe
		
		#now for intra-array we pick PMT=700 to calcuate the intra-array variance, 
		#since PMT=600 is recommended by the manufacturer.
		#cn<-colnames(sboner_cDtf)
		sboner_pmtHiIndex<-col.indices
		
		#go through all these arrays to calculate intra-arra variances and pool together
		#go through the first block to get all the indices of proteins
		sboner_blk<-unique(sboner_cgDtf[,"Block"])
		#rw<-unique(cgDtf[cgDtf[,"Block"]==1,"Row"])
		#cl<-unique(cgDtf[cgDtf[,"Block"]==1,"Column"])
		sboner_cblkOne<-sboner_cgDtf[sboner_cgDtf[,"Block"]==1,]
		#eblkOne<-egDtf[egDtf[,"Block"]==1,]
		
		sboner_ccvs<-data.frame()
		#ecvs<-data.frame()
		
		#do summarizing and plotting
		for(i in sboner_pmtHiIndex)
		{
			#a<-cDtf[,i]
			#go through blocks 1
			sboner_ccv_dtf<-data.frame(cv=rep(0,dim(sboner_cblkOne)[1]),array=rep(i,dim(sboner_cblkOne)[1]))
			sboner_ccv_dtf<-cbind(sboner_ccv_dtf,sboner_cblkOne[,c("Row","Column","Description","Name")])
			
			##summary the data, control first
			for(j in 1:length(sboner_cblkOne[,1]))
			{
				crw<-sboner_cblkOne[j,"Row"]
				ccl<-sboner_cblkOne[j,"Column"]
				#then look through to get CV
				centry<-sboner_cDtf[sboner_cgDtf[,"Row"]==crw&sboner_cgDtf[,"Column"]==ccl,i]
				cgenes<-sboner_cgDtf[sboner_cgDtf[,"Row"]==crw&sboner_cgDtf[,"Column"]==ccl,c("Description","Name")]
				if(length(unique(cgenes)[,1])!=1)
				{
					stop("the genes for variance calculatoin are not grouped correctly.")
				}
				ccv<-abs(sqrt(var(centry))/mean(centry))
				sboner_ccv_dtf[j,"cv"]<-ccv
				sboner_ccv_dtf[j,"Row"]<-crw
				sboner_ccv_dtf[j,"Column"]<-ccl
				sboner_ccv_dtf[j,"Name"]<-cgenes[1,"Name"]
				sboner_ccv_dtf[j,"Description"]<-cgenes[1,"Description"]
				#cv_dtf
				
				#cv_dtf
			}
			if(length(sboner_ccvs)==0)
			{
				sboner_ccvs<-sboner_ccv_dtf
			}else{
				sboner_ccvs<-rbind(sboner_ccvs,sboner_ccv_dtf)
			}

		}#end of for loop
		
		#now return 
		sboner_ccvs
	}
	
	#calculate intra-array variance and mean
	#input:
	#	datamatrix, datamatrix holding the signals for proteins/genes
	#			rows are genes and columns are arrays
	#			also assuming the matrix has been correctly aggregated
	#			most likely this is the control protein fluorescence matrix
	#	geneframe, data frame holding the genes information
	#			it is also aggregated together with datamatrix
	#	col.indices, this is the indices for arrays to be included in
	#			the analysis. If missing, assuming all columns to be used
	#### '   @return the data frame holding the CVs for the control 
	#'@export
	intraArrayMVs<-function(datamatrix, geneframe, col.indices
	)
	{
		if(missing(datamatrix))
		{
			stop("please specify the input data frame.")
		}
		if(missing(geneframe))
		{
			stop("please specify the input gene information frame.")
		}
		if(missing(col.indices))
		{#means doing all
			col.indices<-c(1:length(dataframe[1,]))
		}
		#start doing the variance intra-array with controls
		#intra-array variance first
		sboner_cDtf<-datamatrix
		sboner_cgDtf<-geneframe
		
		#now for intra-array we pick PMT=700 to calcuate the intra-array variance, 
		#since PMT=600 is recommended by the manufacturer.
		#cn<-colnames(sboner_cDtf)
		sboner_pmtHiIndex<-col.indices
		
		#go through all these arrays to calculate intra-arra variances and pool together
		#go through the first block to get all the indices of proteins
		sboner_blk<-unique(sboner_cgDtf[,"Block"])
		#rw<-unique(cgDtf[cgDtf[,"Block"]==1,"Row"])
		#cl<-unique(cgDtf[cgDtf[,"Block"]==1,"Column"])
		sboner_cblkOne<-sboner_cgDtf[sboner_cgDtf[,"Block"]==1,]
		#eblkOne<-egDtf[egDtf[,"Block"]==1,]
		
		sboner_ccvs<-data.frame()
		#ecvs<-data.frame()
		
		#do summarizing and plotting
		for(i in sboner_pmtHiIndex)
		{
			#a<-cDtf[,i]
			#go through blocks 1
			sboner_ccv_dtf<-data.frame(var=rep(0,dim(sboner_cblkOne)[1]),mean=rep(0,dim(sboner_cblkOne)[1]),array=rep(i,dim(sboner_cblkOne)[1]))
			sboner_ccv_dtf<-cbind(sboner_ccv_dtf,sboner_cblkOne[,c("Row","Column","Description","Name")])
			
			##summary the data, control first
			for(j in 1:length(sboner_cblkOne[,1]))
			{
				crw<-sboner_cblkOne[j,"Row"]
				ccl<-sboner_cblkOne[j,"Column"]
				#then look through to get CV
				centry<-sboner_cDtf[sboner_cgDtf[,"Row"]==crw&sboner_cgDtf[,"Column"]==ccl,i]
				cgenes<-sboner_cgDtf[sboner_cgDtf[,"Row"]==crw&sboner_cgDtf[,"Column"]==ccl,c("Description","Name")]
				if(length(unique(cgenes)[,1])!=1)
				{
					stop("the genes for variance calculatoin are not grouped correctly.")
				}
				vars<-var(centry);means<-mean(centry)
				sboner_ccv_dtf[j,"var"]<-vars
				sboner_ccv_dtf[j,"mean"]<-means
				sboner_ccv_dtf[j,"Row"]<-crw
				sboner_ccv_dtf[j,"Column"]<-ccl
				sboner_ccv_dtf[j,"Name"]<-cgenes[1,"Name"]
				sboner_ccv_dtf[j,"Description"]<-cgenes[1,"Description"]
				#cv_dtf
				
				#cv_dtf
			}
			if(length(sboner_ccvs)==0)
			{
				sboner_ccvs<-sboner_ccv_dtf
			}else{
				sboner_ccvs<-rbind(sboner_ccvs,sboner_ccv_dtf)
			}

		}#end of for loop
		
		#now return 
		sboner_ccvs
	}
	
	
	#'@export
	interArrayCVs<-function(cmatrix, cgenes, col.indices)
	{
		
		if(missing(col.indices))
		{#means doing all
			col.indices<-c(1:length(cmatrix[1,]))
		}
		if(missing(cmatrix))
		{
			stop("please specify the input data frame.")
		}
		if(missing(cgenes))
		{
			stop("please specify the input gene information frame.")
		}
		cv<-abs(sqrt(apply(cmatrix[,col.indices], 1, var))/apply(cmatrix[,col.indices], 1, mean))
		dtf<-cbind(cv,cgenes);

		dtf
	}
	
		#'@export
	interArrayMVs<-function(cmatrix, cgenes, col.indices)
	{
		
		if(missing(col.indices))
		{#means doing all
			col.indices<-c(1:length(cmatrix[1,]))
		}
		if(missing(cmatrix))
		{
			stop("please specify the input data frame.")
		}
		if(missing(cgenes))
		{
			stop("please specify the input gene information frame.")
		}
		cvar<-apply(cmatrix[,col.indices], 1, var)
		cmean<-apply(cmatrix[,col.indices], 1, mean)
		dtf<-cbind(cvar,cmean,cgenes);
		colnames(dtf)[1:2]<-c("variance","mean")
		dtf
	}
	
	###function to get rid of small flot data set
	##updated 11/26/2017 to add compare all 
	##input:
	## mtx, data matrix, containing the repeated data in rows next to each other (non-aggregated)
	##		or aggregated data. The rows are for gene entries and columns for different array
	## aggregated, is indicating whethere the data have been aggregated or not. In case of the latter,
	##		the duplicated data are one next to each other. For current implementation, we remove the both 
	##		duplicated entries when either one has smaller signal strength. Can do other ways in the future.
	## threshold, below which the data will be removed
	## index, used to indicate for which row (single row), we want to compare the threshold with
	##		if specified, we usually specify the column with possibly the larget signal, assuming
	##		all other columns will have smaller values. If not specified, we will compare 
	##		signals from all arrays(columns) and pick the gene entries to remove where all 
	##		arrays have smaller signal than thresholds.
	## 		If Index is specified in this case, it will override the method input
	## method, used to define the way to compare the threshold when index is not set. 
	#			ALL, all values has been lower than threshold
	###			Any, if any number in the data row smaller than threshold, the whole row will be removed
	###			Mean, take the mean of the whole row and compare with the thresheld
	##			Median, take median and compare.
	##return:
	##		list contain data matrix cleaned of low signal entries and indice removed according to original matrix
	#'@export
	rmLows<-function(mtx, threshold=40, aggregated=F, index, method=c("ALL", "ANY", "MEAN", "MEDIAN"))
	{
		if(missing(mtx))
		{
			stop("ERROR:please specify the input data matrix")
		}
		if(class(mtx)!="matrix")
		{
			stop("ERROR:expecting data matrix only.")
		}
		flag_compare_all<-F
		if(missing(index))
		{
			flag_compare_all<-T
		}
		
		idx<-c()
		value_by_row<-mtx[,1]
		if(flag_compare_all)
		{
			method<-match.arg(method)
			
			switch(method,
				ALL={value_by_row<-apply(mtx,1, max)},
				MEAN={value_by_row<-apply(mtx,1, mean)},
				MEDIAN={value_by_row<-apply(mtx,1, median)},
				ANY={value_by_row<-apply(mtx,1, min)}
			)
			#now we need to figure out the indice to remove
			idx<-which(value_by_row<threshold)
			#idx_logic<-mtx<threshold
			#for(i in 1:length(mtx[,1]))
			#{
			#	if(sum(idx_logic[i,])==length(idx_logic[i,]))
			#	{
			#		idx<-c(idx,i)
			#	}
			#}
		} else
		{
			idx<-which(mtx[,index]<threshold)
			
		}	
		if(!aggregated)  #need to figure out the indices
		{
			idx<-unique(ceiling(idx/2))
			idx<-rep(idx,rep(2,length(idx)))
			idx<-idx*2
			dif<-c(1,0)
			dif<-rep(dif, length(idx)/2)
			idx<-idx-dif
		}
		rmtx<-mtx
		if(length(idx)>0)
		{
			rmtx<-mtx[-(idx),]
		}
		#else{
		#	mtx
		#}
		
		list(m=rmtx, indices=idx)
	}
	#'@title remove low signal features
	#'@description remove the features whose signal level lower than 
	#'	certain preset threshold from control.
	#'
	#'@details this is the function to rm entries by genes, meaning
	#' we need to somehow go look at the genes entries across
	#' all the arrays and delete the genes entries from
	#' all arrays. This is very specific for control proteins.
	#' It check all the genes across the arrays and all the blocks.
	#' If satisfy the criteria, delete the genes entries from all arrays and all blocks.
	#'
	#' Currently, the criteria is simply as all the means of the gene by each array are 
	#' smaller than the threshold.
	#'	It doesn't make any difference whether the genes are aggregated or not.
	#' Further, we don't implement other methods.
	#'	It intends to do it on control protein field, but could be used for 
	#'	feature field, E. The latter use needs testing.
	#'@param cmtx data matrix for the control protein field in the original 
	#'	elist read from protoArray, eg elist.protoArray$C
	#'@param genes data frame table containing the control gene meta data, elist.protoArray$cgenes
	#'@param threshold numeric the threshold lower than which the control gene entries will
	#'	be removed from the matrix. It is in log scale and 6 by default 
	#'
	#'@return data list contains 3 fields, C, the matrix contains the entries 
	#'	whose signal levels are greater than the threshold; cgene, the relative 
	#'	cgene meta data frame; indices, the indices of entries in the original
	#'	input control protein table that have been removed.
	#
	#'@export
	rmLowCtlGenes<-function(cmtx, genes, threshold=6 
			#,aggregated=T#, method=c("ALL", "MEDIAN", "MEAN", "ANY")
		)
	{
		if(missing(cmtx))
		{
			stop("ERROR:please specify the input data matrix")
		}
		if(class(cmtx)!="matrix")
		{
			stop("ERROR:expecting data matrix only.")
		}
		
		if(missing(genes))
		{
			stop("ERROR:please specify the input gene data frame")
		}
		if(class(genes)!="data.frame")
		{
			stop("ERROR:expecting data frame for gene table only.")
		}
		
		if(dim(cmtx)[1]!=dim(genes)[1])
		{
			stop("ERROR: data matrix and data frame gene table don't have equal number of entries")
		}
		
		#method<-match.arg(method)
		#done and ready for analysis
		
		#first, we need to check the gene table to get all the genes
		unique_genes<-unique(genes[,"Name"])
		
		row_index_to_remove<-c()
		for(i in 1:length(unique_genes))
		{
			#first get all the signal entries for these genes
			index_genes<-which(genes[,"Name"]==unique_genes[i])
			
			#array by gene
			mtx_by_gene<-cmtx[index_genes,]
			mtx_by_gene_mean<-apply(mtx_by_gene,2,mean)
			if(sum(mtx_by_gene_mean<threshold)==length(mtx_by_gene_mean))
			{
				row_index_to_remove<-c(row_index_to_remove,index_genes)
			}
		}
		
		#if we are here, we are done
		
		rmtx<-cmtx;
		rgenes<-genes;
		if(length(row_index_to_remove)>0)
		{
			row_index_to_remove<-sort(row_index_to_remove)
			rmtx<-cmtx[-(row_index_to_remove),]
			rgenes<-genes[-(row_index_to_remove),]
		}
		#else{
		#	mtx
		#}
		
		list(C=rmtx, cgenes=rgenes, indices=row_index_to_remove)
	}
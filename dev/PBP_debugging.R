####this is the code to fix the low expression data
##using anti-histone and isotype control data

#---updated 3/15/2019 Feng@BU

#++++++++++++++

#load libraries
library(PAA)
library(PAST)
library(limma)
library(minpack.lm)

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
	
	ad.elistBC.hf<-ad.elistBC 
	ad.elistBC.hf$E<-ad.elistBC$E[,c(19:36,1:18)]
	ad.elistBC.hf$C<-ad.elistBC$C[,c(19:36,1:18)]
	ad.elistBC.hf$targets<-ad.elistBC$targets[c(19:36,1:18),]

	x<-pmts
	object<-ad.elistBC.hf
#now get one array to do testing.
array.num<-0
index.firstArr<-c(1:6)+6*array.num;
object$E<-object$E[,index.firstArr];
object$C<-object$C[,index.firstArr];
object$targets<-object$targets[index.firstArr,];
object.oneArray<-object
#now getting a blocked data to do fitting.
num.block<-12
block.start<-37
cgenes<-object$cgenes;
genes<-object$genes;

index.block<-which(is.element(cgenes[,"Block"],c(block.start:(block.start+num.block-1))))
object$C<-object$C[index.block,]	
object$cgenes<-object$cgenes[index.block,]

index.block<-which(is.element(genes[,"Block"],c(block.start:(block.start+num.block-1))))
object$E<-object$E[index.block,]
object$genes<-object$genes[index.block,]

B<-log(1.1)
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
inx<-145;

#inx<-124;
y<-object$C
y<-gainAdjust.aggregate(y, "geometric")
#y<-rbind(y,y);
fit_phi<-gainAdjust_fit_Pbp_phi(object=y,x=x, B=ftf$par[1], delta=ftf$par[2], debug=T)

#make a list of para
count<-length(fit_phi);
fit_par<-list()
fit_par[["par"]]<-c(ftf$par[1:2],fit_phi[order(fit_phi)[1:count]])
gainAdjust_plot_Pbp_raw(LMfit=fit_par,y=y[order(fit_phi)[1:count],], x=x, data.aggregated=T,
			G=pmt_saturated,log="xy"
			)
gainAdjust_plot_Pbp(LMfit=fit_par,y=y[order(fit_phi)[1:count],], x=x, data.aggregated=T,
			G=pmt_saturated, log="xy",type="l"
			)
			
#now start plot the data for the high prob/sigal to noise ratio ones and then see what make them up and down
array.num<-c(1:6);

plot(c(350, 700),c(0.1, 2^16-1),log="xy",type="n")
gene.inx<-9540;#8014;#7315;#5618;#3982;#7239;#8798;#8276;
for(i in array.num)
{
	lines(x,object$E[gene.inx*2-1,(i-1)*6+c(1:6)], type="b", col=ceiling(i), pch=ceiling(i/3))
	lines(x,object$E[gene.inx*2,(i-1)*6+c(1:6)], type="b", col=ceiling(i), pch=ceiling(i/3))
	
}


s.norm$E[gene.inx]
s.raw$E[gene.inx]
tts.norm$E$statistic[gene.inx]
tts.norm$E$p.value[gene.inx]
tts.raw$E$statistic[gene.inx]
tts.raw$E$p.value[gene.inx]
tts.pbp.norm$E$statistic[gene.inx]
tts.pbp.norm$E$p.value[gene.inx]

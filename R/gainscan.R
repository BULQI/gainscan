##################################################
# A R Software Tool for Protein Array analysis
#    -1.use to read/import the data
#	 -2.background correct the raw data
#    -3.normalize between arrays, by a linear model on log-transformed data
#    -4.statistically identify the differentially bound proteins between arrays
#    	corrected for multiple comparison
#	 -5.test the shifted distribution/pattern of Ab binding in comparison with the control array
#            Kolmogorov-Smirnov Test or Chi-square test
#	adopted from ARPPA package
#      developed by Feng Feng @ Boston University. 
#      All right reserved
#          10/27/2017
####################################################

##the following is the code to install necessary libs manually
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
#####################

##########Dependency###########
##  library(limma)
####################

#' @include fileIO.R
# # ' @include data.R

#'
#' @title gainscan allows to preprocess and analyze protoarray data
#'
#' @description An R package to read, preprocess and analyze protoarray data
#'
#' @details The functions you're likely to need from \pkg{gainscan} are
#' to read and normalize arrays. Most importantly we can process the data
#' read multiple times under different PMT voltages. The data from multiple
#' scans will be fit for a power function with a baseline and then the incident
#' light signal from each spot will be estimated. This way, the data quality is
#' improved.\cr
#' Otherwise refer to the vignettes to see
#' how to format the documentation.
#'\cr
#'-------------\cr
#'ToDo:
#'\cr
#'1)To restructor the p_elist fitting return object, adding the part for the control object. 
#'	because it could happen that the control have different parameters, when fitted differently. 3/24/2019
#'\cr
#'2)need to see whether to we need to "acquire" the data differently. meaning using "normalized"
#'	parameters.????3/24/2019
#'\cr Answer:No!!!! We are not going to do this. We always do fitting and then normalization RLM, 
#'	but not in an opsite order, RLM-> multi-gain fitting(????) 
#'\cr
#'3)need to check to see array.aggregate vs. gainAdjust.aggregate functions 3/19/2019
#'\cr Answer: they almost doing the same thing, but array.aggregate is a little better doing more?? (4/8/2019)
#'\cr
#'4)need to test/debug data loading/reading modules.v
#'\cr Answer: done with debugging. seems working. 4/9/2019


"_PACKAGE"


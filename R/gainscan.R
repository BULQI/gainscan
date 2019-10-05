##################################################
# An R Software Tool for Protein Microarray Data Analysis
#    -1.use to read/import the data
#	 -2.background correct the raw data
#	 -3.apply the gainscan power funciton model to integrate the array data acquired under different PMT settings
#    -3.normalize between arrays, by a linear model on log-transformed data
#  ToDo:
#    -4.statistically identify the differentially interacting protein features between samples
#    	with correction for multiple comparison
#	 
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
#' @title gainscan to preprocess and analyze protoarray data
#'
#' @description An R package to read, preprocess and analyze protoarray data
#'
#' @details The functions you are likely to need from \pkg{gainscan} are
#' to read and normalize arrays. Most importantly we can process array data
#' acquired multiple times under different PMT voltages. The data from multiple
#' scans will be fitted for a power function with a baseline and then the incident
#' light signal from each spot will be estimated. This way, the data quality is
#' improved and the technical variations are significantly reduced.\cr
#' Please refer to the vignettes to see the explanations and examples.
#' 
#'\cr
#'-------------\cr
#'Currently the package can do the followings \cr
#' \enumerate{
#' \item use to read/import the data
#' \item background correct the raw data
#' \item apply the gainscan power funciton model to integrate the array data acquired under different PMT settings
#' \item normalize between arrays, by a linear model on log-transformed data \cr\cr
#'ToDo:
#' \item run differential analysis to statistically identify the positive interacting protein features 
#'    	with correction for multiple comparison
#' \item add a section to describe in vignettes about the differential analysis and FDR.
#'	}
#'\cr
#'


"_PACKAGE"


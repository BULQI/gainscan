## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")

library(qvalue)

###note for a uniform distributed p values, p.adj is identical to qvalues (??why)
p<-runif(1000)

p_sort<-sort(p)

p_BH<-p.adjust(p_sort, method="BH")

p_q<-qvalue(p_sort)

p_BH
p_q$qvalue

########example from R documentation, non-uniform distributed p values
# import data
data(hedenfalk)
p <- hedenfalk$p

p_sort<-sort(p)
# get q-value object
qobj <- qvalue(p_sort)
plot(qobj)
hist(qobj)

qobj$qvalue[c(1:100)]

p_BH<-p.adjust(p_sort, method="BH")
p_BH[c(1:100)]
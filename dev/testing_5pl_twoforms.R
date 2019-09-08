#Test how to make the out of limit case of 5pl to get solved.

#two forms of 5pl, log vs regular (in terms of x)
#	regular: y=a + (d-a)/(1+(b/x)^c)^g
#				a, lower; d, upper; b, mid point; c, slope of mid range and g, asymetric factor

#	log: y=a+(d-a)/(1+exp((xmid-x)/scal))^g

# to transform them, we have
#	y=a+(d-a)/(1+exp(log(b/x)^c))^g
#	y=a+(d-a)/(1+exp(c*log(b)-log(x)))^g
#   y=a+(d-a)/(1+exp((log(b)-log(x))/(1/c)))^g
#		so between them xmid==log(b), scal==1/c

##we are looking for a solution to overcome the issue where we run out of limit
## when calculate exp((xmid-x)/scal)
## the solution is like this,
## take a exp(log((d-a)/(1+exp((xmid-x)/scal))^g))
##   rearrange this  exp(log(d-a)-log((1+exp(xmid-x)/scal)^g))
##					exp(log(d-a)-g*log(1+exp(xmid-x)/scal))
##		there are one case we need to consider, when exp is too big or exp is too small (both relative to 1).
##		where exp() is too big or (xmid-x)>709, we approximate to get log(1+exp(xmid-x)/scal)~=(xmid-x)/scal 
##		the whole thing becomes exp(log(d-a)-g*(xmid-x)/scal)
fpl.regular<-function(a, d, b,c,g,x)
{
	y<-a + (d-a)/(1+(b/x)^c)^g
}


fpl.log<-function(a, d, xmid,scal,g,x)
{
	y<-a+(d-a)/(1+exp((xmid-x)/scal))^g
}
fpl.log.app<-function(a, d, xmid,scal,g,x)
{
	y<-x; #make the holder first
	y<- -1
	#split the array into to parts
	regular<-which((xmid-x)/scal<=709)
	y[regular]<-a+(d-a)/(1+exp((xmid-x[regular])/scal))^g
	y[-regular]<-a+exp(log(d-a)-g*(xmid-x[-regular])/scal)	
	y
}

x<-seq(4.5,8,by=0.01)
a<-0.1
d<-65555
xmid<-6.5
scal<-0.002
g<-0.000807

y<-fpl.log(a,d,xmid,scal,g,x)
plot(x,y, type="l")

y<-fpl.regular(a,d, exp(xmid),1/scal,g,exp(x))
lines(x,y,col=2,lty=2)

y<-fpl.log.app(a,d, xmid,scal,g,x)
lines(x,y,col=3,lty=4, lwd=2)
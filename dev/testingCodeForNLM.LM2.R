library(minpack.lm)
> ??mCall
> SSfpl
function (input, A, B, xmid, scal) 
{
    .expr1 <- B - A
    .expr2 <- xmid - input
    .expr4 <- exp(.e2 <- .expr2/scal)
    .expr5 <- 1 + .expr4
    .value <- A + .expr1/.expr5
    .actualArgs <- as.list(match.call()[c("A", "B", "xmid", "scal")])
    if (all(unlist(lapply(.actualArgs, is.name)))) {
        .expr8 <- 1/.expr5
        .expr13 <- .expr5^2
        .grad <- array(0, c(length(.value), 4L), list(NULL, c("A", 
            "B", "xmid", "scal")))
        .grad[, "A"] <- 1 - .expr8
        .grad[, "B"] <- .expr8
        .grad[, "xmid"] <- -(xm <- .expr1 * .expr4/scal/.expr13)
        .grad[, "scal"] <- xm * .e2
        dimnames(.grad) <- list(NULL, .actualArgs)
        attr(.value, "gradient") <- .grad
    }
    .valueNLS
}
<bytecode: 0x000000000c6f1118>
<environment: namespace:stats>
attr(,"initial")
function (mCall, data, LHS) 
{
    xy <- sortedXyData(mCall[["input"]], LHS, data)
    if (nrow(xy) < 5) {
        stop("too few distinct input values to fit a four-parameter logistic")
    }
    rng <- range(xy$y)
    drng <- diff(rng)
    xy$prop <- (xy$y - rng[1L] + 0.05 * drng)/(1.1 * drng)
    ir <- as.vector(coef(lm(x ~ I(log(prop/(1 - prop))), data = xy)))
    pars <- as.vector(coef(nls(y ~ cbind(1, 1/(1 + exp((xmid - 
        x)/exp(lscal)))), data = xy, start = list(xmid = ir[1L], 
        lscal = log(abs(ir[2L]))), algorithm = "plinear")))
    value <- c(pars[3L], pars[3L] + pars[4L], pars[1L], exp(pars[2L]))
    names(value) <- mCall[c("A", "B", "xmid", "scal")]
    value
}
<environment: namespace:stats>
attr(,"pnames")
[1] "A"    "B"    "xmid" "scal"
attr(,"class")
[1] "selfStart"
> .expr1
Error: object '.expr1' not found



## weighted nonlinear regression
Treated <- Puromycin[Puromycin$state == "treated", ]
conc<-Treated$conc
Vm<-6
K<-10
resp<-(Vm*conc)/(K+conc)

simData<-data.frame("conc"=conc, "resp"=resp)
weighted.MM <- function(resp, conc, Vm, K)
{
    ## Purpose: exactly as white book p. 451 -- RHS for nls()
    ##  Weighted version of Michaelis-Menten model
    ## ----------------------------------------------------------
    ## Arguments: 'y', 'x' and the two parameters (see book)
    ## ----------------------------------------------------------
    ## Author: Martin Maechler, Date: 23 Mar 2001

    pred <- (Vm * conc)/(K + conc)
    (resp - pred) #/ sqrt(pred)
}

##

Pur.wt <- nlsLM( ~ weighted.MM(resp, conc, Vm, K), data = simData,
              start = list(Vm = 1, K = 1))
summary(Pur.wt)

###########----------------
## values over which to simulate data
x <- seq(0,5,length=100)
## model based on a list of parameters
getPred <- function(parS, xx) parS$a * exp(xx * parS$b) + parS$c
## parameter values used to simulate data
pp <- list(a=9,b=-1, c=6,k=6)
## simulated data, with noise
simDNoisy <- getPred(pp,x) + rnorm(length(x),sd=.1)
## plot data
plot(x,simDNoisy, main="data")
## residual function
residFun <- function(p, observed, xx) observed - getPred(p,xx)
## starting values for parameters
parStart <- list(a=3,b=-.001, c=1)
## perform fit
nls.out <- nls.lm(par=parStart, fn = residFun, observed = simDNoisy,
xx = x, control = nls.lm.control(nprint=1))
## plot model evaluated at final parameter estimates
lines(x,getPred(as.list(coef(nls.out)), x), col=2, lwd=2)

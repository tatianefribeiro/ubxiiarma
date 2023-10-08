# Implemented by Fabio M Bayer (bayer@ufsm.br) em 15/10/2015
# Changes on May 17, 2015
#
# Some informations:
# diag = 0 : It does not plot the graphics 
# diag = 1 : It plots the graphics on the screen
# diag = 2 : It generates the plots in PDF and plot the graph on the screen
#
# h : number of steps ahead to make a forecast
#
# The input object of the function should be a time series (ts)
#
# There are four types of residuals to be used, with 'resid' ranging from 1 to 4.
#
# Usage examples:
#
# 1) BARMA(2,3) with logit link function and h = 6 
# fit <- barma(y,ar=c(1,2),ma=c(1,2),resid=2)
# 
# Obs: Note that you can choose any lags you desire.
# 
# 2) Printing graphs in PDF files with probit link function
# fit <- barma(y,ar=c(1,2),ma=c(1,2),diag=2,link="probit",resid=3)


barma<- function (y, ar=NA, ma=NA, link = "logit",diag=1,h=6,X=NA,X_hat=NA,resid=1)
{  
  source("barma.fit.r")
  
  if (min(y) <= 0 || max(y) >= 1)
    stop("OUT OF RANGE (0,1)!")
  
  if(is.ts(y)==T)
  {
    freq<-frequency(y)
  }else stop("data can be a time-series object")
  
  
  if(any(is.na(ar))==F) names_phi<-c(paste("phi",ar,sep=""))
  
  if(any(is.na(ma))==F) names_theta<-c(paste("theta",ma,sep=""))
  
  if(any(is.na(X))==F)
  {
    names_beta<-c(paste("beta",1:ncol( as.matrix(X) ),sep=""))
  }
  
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  
  m <- max(p,q,na.rm=T)
  
  p1 <- length(ar)
  q1 <- length(ma)
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))
    stats <- make.link(linktemp)
  else stop(paste(linktemp, "link not available, available links are \"logit\", ",
                  "\"probit\" and \"cloglog\""))
  
  link1 <- structure(list(link = linktemp, 
                          linkfun = stats$linkfun,
                          linkinv = stats$linkinv, 
                          mu.eta = stats$mu.eta, 
                          diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t)))
  )
  )
  
  
  fit1 <- barma.fit(y, ar, ma, link1, names_phi, names_theta, names_beta, diag, h, X, X_hat,resid=resid) # model estimation
  
  return(fit1)
}


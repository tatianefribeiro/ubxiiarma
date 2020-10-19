simu.UBXIIarma <- function(n,phi=NA,theta=NA, alpha=1,c_par=2.2,tau=0.5,freq=12,link="logit")
{
  source("UBXII-funcs.R")
  if(any(is.na(phi)==F))
  {
    ar <- 1:length(phi)
  }else{
    ar<-0
    phi<-0
  }
  
  if(any(is.na(theta)==F))
  {
    ma <- 1:length(theta)
  }else{
    ma<-0
    theta<-0
  }
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))
  {  
    stats <- make.link(linktemp)
  }else{
    stop(paste(linktemp, "link not available, available links are \"logit\", ",
               "\"probit\" and \"cloglog\""))
  } 
  
  link <- structure(list(link = linktemp, 
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv
  )
  )
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  
  {
    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)  ###!!!!!?? pg.390, KARMA?
    
    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)
    
    error<-rep(0,n+m) 
    eta<- y <- rep(NA,n+m)  # PQ ATÉ n+m?
    
    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + as.numeric(phi%*%ynew[i-ar]) + as.numeric(theta%*%error[i-ma])
      mu[i]   <- linkinv(eta[i])
      y[i]    <- r_UBXII(1,mu[i],c_par)
      ynew[i] <- linkfun(y[i])       #!!Onde considero o componente aleotório?? (Rayleigh nesse caso.)
      error[i]<- ynew[i]-eta[i]   
    }
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )  #pq até n+m e não até m, como é feito no .fit
  } 
}


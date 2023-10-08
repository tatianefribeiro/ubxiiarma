# Created by Fernando A. Peña-Ramírez, Renata R. Guerra and Tatiane F. Ribeiro (tatianefr@ime.usp.br), October/2020
# Based on  barma.fit.r codes created by Fabio M. Bayer on October 15, 2015
#
# Some informations:
# diag = 0 : It does not plot the graphics 
# diag = 1 : It plots the graphics on the screen
# diag = 2 : It generates the plots in PDF and plot the graph on the screen
#
# h : number of steps ahead to make a forecast

ubxiiarma.fit<-function (y, ar = NA, ma = NA, tau = .5,link = "logit",
                         h=0, diag=0,X = NA,X_hat=NA)
{
  source("ubxii-funcs.r")
  if (min(y) <= 0 || max(y) >= 1)
    stop("OUT OF RANGE (0,1)!")
  
  if(is.ts(y)==T)  freq<-frequency(y) else stop("data can be a time-series object")
  
  z<-c()
  maxit1<-50
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  y1 <- y[(m+1):n]
  p1 <- length(ar)  
  q1 <- length(ma)
  error <- rep(0,n) 
  eta <- rep(NA,n)
  
  y_prev <- c(rep(NA,(n+h))) 
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog"))){
    stats <- make.link(linktemp)
  }  else {
    stop(paste(linktemp, "link not available, available links are \"logit\", ","\"probit\" and \"cloglog\""))
  }
  
  link = linktemp 
  linkfun = stats$linkfun
  linkinv = stats$linkinv 
  mu.eta = stats$mu.eta 
  diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t)))
  ynew = linkfun(y) 
  ynew_ar <- suppressWarnings(matrix(ynew,(n-1),max(p,1,na.rm=T)))
  
  ###########################################################
  if(any(is.na(ar)) == F) {
    names_phi <- c(paste("phi", ar, sep = ""))
    Z <- suppressWarnings(matrix(ynew, (n-1), p1)[m:(n-1),])} else {
      ar = p1<-0; Z <- NA  
    } 
  
  if(any(is.na(ma)) == F) {
    names_theta <- c(paste("theta", ma, sep = ""))
  } else ma = q1 <- 0 
  
  if(any(is.na(X)) == F){
    names_beta<-c(paste("beta", 1 : ncol(as.matrix(X)), sep = ""))
    Xm <- X[(m+1):n, ]     
    k = ncol(X)
  } else {
    k = 0 
    X <- matrix(rep(0,n), nrow = n)
    Xm <- NA
  }
  
  # Recorrences
  q_1 <- max(q1, 1)
  R <- matrix(rep(NA, (n-m)*q_1), ncol = q_1)  
  k_i <- q1/q_1                                 
  
  deta.dalpha <- rep(0, n)
  deta.dbeta <- matrix(0, ncol=max(k,1), nrow=n)
  deta.dphi <- matrix(0, ncol=p1, nrow=n)
  deta.dtheta <- matrix(0, ncol=q_1, nrow=n)
  
  Xstart <- (cbind(rep(1, (n-m)), Xm, Z))         
  Xstart <- matrix(apply(Xstart, 1, na.omit),nrow = (n-m),byrow = T)
  ols <- lm.fit(Xstart, ynew[(m+1) : n])$coef
  initial <- c(rep(0, k+p1+q1+1),1)
  initial[1 : (k+p1+1)] <- ols
  
  loglik <- function(z) 
  {
    alpha <- z[1]
    if(k==0)  beta = as.matrix(0) else beta = as.matrix(z[2:(k+1)])
    if(p1==0) {phi = as.matrix(0);ar=1} else phi = as.matrix(z[(k+2):(k+p1+1)]) 
    if(q1==0) theta = as.matrix(0) else  theta = as.matrix(z[(k+p1+2):(k+p1+q1+1)])
    c_par <- z[length(z)]
    
    Xbeta <- X%*%beta
    Xbeta_ar <- suppressWarnings(matrix(Xbeta, (n-1), max(p, 1, na.rm = T)))
    
    for(i in (m+1):n)
    {
      eta[i] <- alpha + Xbeta[i] + (ynew_ar[(i-1), ar] - Xbeta_ar[(i-1), ar])%*%phi + t(theta)%*%error[i-ma]
      error[i] <- ynew[i] - eta[i] 
    }
    q_t <- linkinv(eta[(m+1):n])
    
    ll <- log(-c_par*log(tau))-log(y1)+(c_par-1)*log(log(1/y1))-
      log(1+(log(1/y1))^c_par)-log(log(1+(log(1/q_t))^c_par))-
      (log(1/tau)*log(1+(log(1/y1))^c_par))/
      log(1+(log(1/q_t))^c_par)
    sum(ll)
  } 
  
  #######################################################################
  escore.UBXIIarma <- function(z)
  {
    alpha <- z[1]
    if(k==0)  beta = as.matrix(0) else beta = as.matrix(z[2:(k+1)])
    if(p1==0) {phi = as.matrix(0);ar=1} else phi = as.matrix(z[(k+2):(k+p1+1)]) 
    if(q1==0) theta = as.matrix(0) else  theta = as.matrix(z[(k+p1+2):(k+p1+q1+1)])
    c_par <- z[length(z)]
    
    Xbeta <- X%*%beta
    Xbeta_ar <- suppressWarnings(matrix(Xbeta, (n-1), max(p, 1, na.rm = T)))
    for(i in (m+1):n)
    {
      eta[i] <- alpha + Xbeta[i] + (ynew_ar[(i-1),ar] - Xbeta_ar[(i-1),ar])%*%phi + t(theta)%*%error[i-ma]
      error[i] <- ynew[i] - eta[i] 
    }
    
    q_t <- linkinv(eta[(m+1):n])
    
    Xbeta <- X%*%beta
    for(i in 1:(n-m)){
      R[i,] <- error[i+m-ma]*k_i}
    
    for(i in (m+1):n)
    {
      deta.dalpha[i] <- 1 - deta.dalpha[i-ma]%*%theta
      deta.dbeta[i,] <- X[i,] - t(phi)%*%X[i-ar,] - t(theta)%*%deta.dbeta[i-ma,]
      deta.dphi[i,] <- ynew_ar[i-ar]- Xbeta[i-ar] - t(theta)%*%deta.dphi[i-ma,]
      deta.dtheta[i,] <- R[(i-m),] - t(theta)%*%deta.dtheta[i-ma,]
    }
    
    v <- deta.dalpha[(m+1):n]
    rM <- deta.dbeta[(m+1):n,]
    rP <- deta.dphi[(m+1):n,]
    rR <- deta.dtheta[(m+1):n,]
    
    mT <- diag(mu.eta(eta[(m+1):n]))
    
    #ell_q
    q_star <- (log(1/q_t))^(c_par-1)/(q_t*(1+(log(1/q_t))^c_par)*log(1+(log(1/q_t))^c_par))
    q_dag <- log(1/tau)*(log(1/q_t))^(c_par-1)/(q_t*(1+(log(1/q_t))^c_par)*(log(1+(log(1/q_t))^c_par))^2)
    y_star <- log(1+(log(1/y1))^c_par)
    
    a_t <- c_par*(q_star-q_dag*y_star)
    
    #ell_c_par
    y_sust <- as.vector(
      1/c_par+log(log(1/y1))-(log(log(1/q_t))*((1+(log(1/q_t))^c_par)-1))/((1+(log(1/q_t))^c_par)*log(1+(log(1/q_t))^c_par))-
        (((1+(log(1/y1))^c_par)-1)*log(log(1/y1)))/(1+(log(1/y1))^c_par)-
        (log(1/tau)*log(1+(log(1/q_t))^c_par)*((1+(log(1/y1))^c_par)^(-1))*((1+(log(1/y1))^c_par)-1)*log(log(1/y1)))/((log(1+(log(1/q_t))^c_par))^2)+
        (log(1/tau)*((1+(log(1/q_t))^c_par)-1)*log(log(1/q_t))*log(1+(log(1/y1))^c_par))/((1+(log(1/q_t))^c_par)*(log(1+(log(1/q_t))^c_par))^2)
    )

    Ualpha <- t(v) %*% mT %*% a_t
    Ubeta <- t(rM) %*% mT %*% a_t
    Uphi <-   t(rP) %*% mT %*% a_t
    Utheta <- t(rR) %*% mT %*% a_t
    Uc <- sum(y_sust)
    
    rval <- c(Ualpha,Ubeta,Uphi,Utheta,Uc)
    return(rval[rval!=0])
  }
  
  #############################################
  
  opt<-optim(initial, loglik, 
             escore.UBXIIarma, 
             method = "BFGS", hessian = TRUE,
             control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
  
  
  if (opt$conv != 0)
  {
    warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
    opt<-optim(initial, loglik, 
               method = "BFGS", hessian = TRUE,
               control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
    }else{
      warning("IT WORKS WITH NUMERICAL GRADIENT!")
    }
  }
  
  z$conv <- opt$conv
  coef <- (opt$par)[1:(p1+q1+k+2)]
  alpha <- coef[1]
  if(k==0) beta=names_beta=NULL else z$beta <- coef[2:(k+1)]
  if(p1==0) phi=names_phi=NULL else z$phi <- coef[(k+2):(k+p1+1)]
  if(q1==0) theta=names_theta=NULL else z$theta <- coef[(k+p1+2):(k+p1+q1+1)]
  z$c_par <- coef[length(coef)]
  
  names_par <- c("alpha",names_beta,names_phi,names_theta,"c_par")
  names(coef)<-names_par
  z$coeff <- coef
  J_inv <- solve(-(opt$hessian))
  z$stderror<-sqrt(diag(J_inv))
  z$zstat <- z$coef/z$stderror
  z$pvalues <- 2*(1 - pnorm(abs(z$zstat)) )
  z$loglik <- opt$value
  z$counts <- as.numeric(opt$counts[1])
  
  # if(any(is.na(X)==F))
  # {
  #   z$k<- (p1+q1+k+2)
  #   z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
  #   z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
  #   z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  # }else{
  #   z$k<- (p1+q1+2)
  #   z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
  #   z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
  #   z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  # }
  
  if(any(is.na(X)==F))
  {
    z$k<- (p1+q1+k+2)
    z$aic <- -2*(z$loglik)+2*(z$k)
    z$bic <- -2*(z$loglik)+log(n)*(z$k)
    z$hq <- -2*(z$loglik)+log(log(n))*(z$k)
  }else{
    z$k<- (p1+q1+2)
    z$aic <- -2*(z$loglik)+2*(z$k)
    z$bic <- -2*(z$loglik)+log(n)*(z$k)
    z$hq <- -2*(z$loglik)+log(log(n))*(z$k)
  }

  model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
  z$model <- model_presentation

 # Predicted values   (NO REGRESSORS)
  if(k==0){                                     
    alpha <- as.numeric(coef[1])
    phi <- as.numeric(coef[2:(p1+1)])
    theta <- as.numeric(coef[(p1+2):(p1+q1+1)])
    c_par <- as.numeric(coef[p1+q1+2])  
    
    z$alpha <- alpha
    z$phi <- phi
    z$theta <- theta
    z$c_par <- c_par
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    if(p1==0) {phi = as.matrix(0);ar=1}
    if(q1==0) {theta = as.matrix(0);ma=1}
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + ynew_ar[(i-1),ar]%*%as.matrix(phi) + (theta%*%errorhat[i-ma])
      errorhat[i]<- ynew[i]-etahat[i] # predictor scale
    }
    
    q_hat <- linkinv(etahat[(m+1):n])   # fitted values 
    
    z$fitted <- ts(c(rep(NA,m),q_hat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    
    # Forecasting
    ynew_prev <- c(ynew,rep(NA,h))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h)
    {
      ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar]) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 
    }
    
    z$forecast <- y_prev[(n+1):(n+h)]
  } else{                              # with REGRESSORS
    X_hat <- as.matrix(X_hat)
    
    alpha <- as.numeric(coef[1])
    beta <- as.numeric(coef[2:(k+1)])
    phi <- as.numeric(coef[(k+2):(k+p1+1)])
    theta <- as.numeric(coef[(k+p1+2):(k+p1+q1+1)])
    c_par <- as.numeric(coef[length(coef)])  
    
    z$alpha <- alpha
    z$beta <- beta
    z$phi <- phi
    z$theta <- theta
    z$c_par <- c_par
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    if(p1==0) {phi = as.matrix(0);ar=1}
    if(q1==0) {theta = as.matrix(0);ma=1}

    Xbeta <- X%*%as.matrix(beta)
    Xbeta_ar <- suppressWarnings(matrix(Xbeta, (n-1), max(p, 1, na.rm = T)))

    for(i in (m+1):n)
    {
      etahat[i] <- alpha + Xbeta[i] + (ynew_ar[(i-1), ar] - Xbeta_ar[(i-1), ar])%*%as.matrix(phi) +
        t(as.matrix(theta))%*%errorhat[i-ma]
      errorhat[i] <- ynew[i] - etahat[i]
    }
    q_hat <- linkinv(etahat[(m+1):n])   # fitted values 
    
    z$fitted <- ts(c(rep(NA,m),q_hat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    
    # Forecasting
    ynew_prev <- c(ynew,rep(NA,h))
    y_prev[1:n] <- z$fitted
    
    X_prev <- rbind(X,X_hat)
    
    for(i in 1:h)
    {
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) +
        (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta))) +
        (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0
    }
    
    z$forecast <- y_prev[(n+1):(n+h)]
    
  }
  
  
  # Quantile residuals 
  z$residuals <- as.vector(qnorm(cdf_UBXII(y[(m+1):n],z$fitted[(m+1):n],z$c_par)))
  residc <- z$residuals 

  # GRAPHICS 
  
  if(diag>0)
  {
    
    print(model_presentation)
    print(" ",quote=F)
    print(c("Log-likelihood:",round(z$loglik,4)),quote=F)
    print(c("Number of iterations in BFGS optim:",z$counts),quote=F)
    print(c("AIC:",round(z$aic,4)," SIC:",round(z$bic,4)," HQ:",round(z$hq,4)),quote=F)
    
    print("Residuals:",quote=F)
    print(summary(residc))
    
    par(mfrow=c(1,1))
    plot(y,type="l",ylab="Serie",xlab="Time",ylim=c(min(y),max(y)))
    lines(z$fitted,col="blue",lty=2)
      legend("topright",c("Observed data","Predicted median"),
             pt.bg="white", lty=c(1,2), bty="n",col=c(1,"blue"))
    
   
    w1<-5
    h1<-4
    
    if(diag>1)
    {
      postscript(file = "resid_v_ind.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1))
        par(mgp=c(1.7, 0.45, 0))
        plot(residc,main=" ",xlab="Index",ylab="Residuals", pch = "+",
             ylim=c(-4,4))
        lines(t,rep(-3,n+h+6),lty=2,col=1)
        lines(t,rep(3,n+h+6),lty=2,col=1)
        lines(t,rep(-2,n+h+6),lty=3,col=1)
        lines(t,rep(2,n+h+6),lty=3,col=1)
      }
      dev.off()

      postscript(file = "adjusted.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1))  
        par(mgp=c(1.7, 0.45, 0))
        plot(y,type="l",ylab="Serie",xlab="Time")
        lines(z$fitted,col=2,lty=2)
        legend("topright",c("Observed data","Predicted median"),
               pt.bg="white", lty=c(1,2), bty="n",col=c(1,2), cex=.8)
        
      }
      dev.off()
      
    }    
  }  
  
  return(z)
}


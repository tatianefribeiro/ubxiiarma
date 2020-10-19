UBXIIARMA.fit<-function (y, ar = NA, ma = NA, X = NA, tau = .5,link = "logit")
{
  if (min(y) <= 0 || max(y) >= 1)
    stop("OUT OF RANGE (0,1)!")
  
  if(is.ts(y)==T)  freq<-frequency(y) else stop("data can be a time-series object")
  
  z<-c()
  maxit1<-10000
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  y1 <- y[(m+1):n]
  p1 <- length(ar)
  q1 <- length(ma)
  error <- rep(0,n) 
  eta <- rep(NA,n)
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
  
###########################################################3  
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

  ###FB  recorrences   #*
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
    #ones <- rep(1,(n-m))
    
    Ualpha <- t(v) %*% mT %*% a_t
    Ubeta <- t(rM) %*% mT %*% a_t
    Uphi <-   t(rP) %*% mT %*% a_t
    Utheta <- t(rR) %*% mT %*% a_t
    Uc <- sum(y_sust)#t(y_sust) %*% ones
    
    rval <- c(Ualpha,Ubeta,Uphi,Utheta,Uc)
    return(rval[rval!=0])
  }
  
  #############################################
  
  opt<-optim(initial, loglik, 
            # escore.UBXIIarma, 
             method = "BFGS", hessian = TRUE,
             control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
  
  
  if (opt$conv != 0)
  {
    warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
    opt<-optim(initial, loglik, 
               #  escore.UBXIIarma, 
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
  z$zstat <- abs(z$coef/z$stderror)
  z$pvalues <- 2*(1 - pnorm(z$zstat) )
  z$loglik <- opt$value
  z$counts <- as.numeric(opt$counts[1])
  
  if(any(is.na(X)==F))
  {
    z$k<- (p1+q1+1+k)
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }else{
    z$k<- (p1+q1+1)
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }
  
  model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
  z$model <- model_presentation
  return(z)
}

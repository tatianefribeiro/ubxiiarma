rm(list = ls())
setwd("/home/tatiane/Insync/tfr1@de.ufpe.br/Google Drive/master_thesis_Tati/5-scripts_cap3/UBXII-ARMA_simu_res")
source("UBXIIarma.fit.R")
source("simu_UBXIIarma.R")

R <- 10
n <- 75
alpha <- 0.5
phi<-c(.6,-.4)
theta <- c(0.4,.1) 
c_par <- 5.6
model <- c(alpha,phi,theta,c_par)

#n= 75; 125; 200; 300

p1<- 1:length(phi)
q1<- 1:length(theta)

#ONLY FOR TEST
#y <- simu.UBXIIarma(n,phi=phi,theta=theta,alpha=alpha,c_par=c_par)

# For the confidence intervals
alpha_erro <- 0.05
quantil <-1-alpha_erro/2
z<-qnorm(quantil)
t<-qt(quantil, n-1)

# To save the results
estim <-ICi<-ICs<- err <- matrix(NA, nrow = R, ncol = length(model))
calpha<-cphi1<-ctheta1<-cphi2<-ctheta2<-cc_par<-0
i<-0

set.seed(5)

### simulation
tempo.inicio = Sys.time()
while(i<R) 
{ 
  y <- simu.UBXIIarma(n,phi=phi,theta=theta,alpha=alpha,c_par=c_par)
  result <- try(UBXIIARMA.fit(y,ar=p1,ma=q1), silent = T)
  if(class(result) == "try-error" ) 
  {
    result$conv<-1
    convergencia<-1
  }
  
  if(result$conv == 0)
  {
    convergencia<-  result$conv
  }else{
    
    warning("FUNCTION DID NOT CONVERGE!")
  }
  
  if(convergencia==0 && sum(is.na(result$stderror))==0)
  {
    i<-i+1
    print(c("i=",i))
    estim[i,] <- result$model[,1]
    err[i,] <- result$model[,2]
    ICi[i,]<-estim[i,]-(z*err[i,])
    ICs[i,]<-estim[i,]+ (z*err[i,])
    
    if (ICi[i,1]<=alpha && ICs[i,1]>=alpha)
    {
      calpha<-calpha+1
    }
    
    if (ICi[i,2]<= phi[1] && ICs[i,2]>=phi[1])
    {
      cphi1<-cphi1+1
    }
    if (ICi[i,3]<= phi[2] && ICs[i,3]>=phi[2])
    {
      cphi2<-cphi2+1
    }
    if (ICi[i,4]<= theta[1] && ICs[i,4]>=theta[1])
    {
      ctheta1<-ctheta1+1
    }
    if (ICi[i,5]<= theta[2] && ICs[i,5]>=theta[2])
    {
      ctheta2<-ctheta2+1
    }
    if (ICi[i,6]<= c_par && ICs[i,6]>=c_par)
    {
      cc_par<-cc_par+1
    }
  } # fim convergencia MC
} #fim loop MC

tempo.fim = Sys.time()
tempo.exec = tempo.fim- tempo.inicio

### mean
m <- apply(estim, 2, mean)
### bias
bias <- (model-m)
### relative percentage bias
biasP <- bias/model *100
### SD
erro <- apply(estim, 2, sd)
### MSE
MSE <- apply(estim, 2, var)+bias^2
### 
CR<-c(calpha,cphi1,cphi2,ctheta1,ctheta2,cc_par)/R
## final results
results <- rbind(m, bias, biasP, erro, MSE,CR)
rownames(results) <- c("Mean", "Bias","RB%", "SE", "MSE","CR")
print(results)
print(tempo.exec)

#save the results
# ttt<-
#   paste0("/home/tatiane/Insync/tfr1@de.ufpe.br/Google Drive/master_thesis_Tati/5-scripts_cap3/UBXII-ARMA/simu_res",n,
#          "_alpha",alpha,
#          "_phi", phi[1],"_",phi[2],
#          "_theta",theta[1],"_",theta[2],
#          "_c",c_par,
#          "_R",R,".RData")
# 
# save.image(ttt)


# Function to test various BARMA models and returns the best one
# according to AIC criterion
#
# Implemented by Fabio M Bayer (bayer@ufsm.br) on October 02, 2014 
#
# Changed by Renata R Guerra (renata.r.guerra@ufsm.br) on November 13, 2018
# Data of input:
# - time series: it is the time series of interest
# - sf:  information about the start and frequency of the time series should follow the following format:
#     ex: sf<-c(start=c(1994,7),frequency=12)  # it informs the start and the frequency of the time series
# - h: number of observations that will be separated from the sample for out-of-sample forecasting
# - pmax, qmax, Pmax, Qmax: maximum orders to be tested. If these values were not informed, the default values below will be used
#
# OBS.: 
#  1- The forecasts for the regressors are obtained from Holt-Winters.
#  2- It uses nested models 

best.karma<-function(serie, sf, h=6, pmax=6, qmax=6, nbest=10,link = "logit",X=NA,X_hat=NA)
{
  source("karma.r")
  #n<-length(serie)
  #out<-seq((n-h+1),n,1)  # para retirar as ultimas h observacoes
  #y<-ts(serie[-out],start=c(sf[1],sf[2]),frequency=sf[3]) # time series without the last observations
  y<-ts(serie,start=c(sf[1],sf[2]),frequency=sf[3])
  
  
  # It initializes the AIC criteria
  fit<-karma(y, ar=1, diag=0,link = link)
  aicmin<-fit$aic
  
  print(aicmin)
  
  model1<-model2<-model3<-model4<-model5<-0
  model6<-model7<-model8<-model8<-model10<-0
  
  best_aic<-rep(Inf,nbest) # It saves the 10 smallest AICs
  #melhores<-rep(0,(nbest)) # It saves the order of the 10 best models
  melhores<-matrix(rep(0,(nbest*3)),ncol=3) # It saves the order of the 10 best models
  colnames(melhores)<-c("p","q","AIC")
  
  tot<-0 # It initializes the count of how many models will be tested.
  bug<-0 # Initializes a counter for how many times the Arima(...) estimation failed
  
  for(p in 0:pmax)
  { 
    for(q in 0:qmax)
    { 
      if(p==0)   ar1<-NA else ar1<-1:p
      if(q==0)   ma1<-NA else ma1<-1:q
      
      if(sum(is.na(c(ar1,ma1)))<2 )
      {   
        print(c(p,q),quote=F)
        fitkarma<-karma(y, ar=ar1, ma=ma1, diag=0,link = link,X=X,X_hat=X_hat)
        tot<-tot+1
        
        if(fitkarma$conv != 0)
        {  
          print(c("NO CONVERGENCE  ",p,q),quote=F)
          bug<-bug+1
          next # Exit this loop and move to the next one
        }          
        if(aicmin>fitkarma$aic) # best model according to AIC
        {  
          aicmin<-fitkarma$aic
          best_model_aic <- fitkarma$model
          print("###########################################")
          print(aicmin)
        }
        if(fitkarma$aic<max(best_aic))
        {
          maximo<-order(best_aic)[nbest]
          best_aic[maximo]<-fitkarma$aic
          melhores[maximo,]<-c(p,q,fitkarma$aic)
          #print(melhores)
        }
        
        
      }
    }
  }
  
  
  print(" ",quote=F)
  print("SELECTED MODEL FROM AIC",quote=F)
  print(best_model_aic,quote=F)
  print(" ",quote=F)
  print(c("Total of tested models =",tot),quote=F)
  print(c("Total of errors in the estimation =",bug),quote=F)
  
  print(" ",quote=F)
  print("THE BEST MODELS",quote=F)
  print(melhores,quote=F)
  
}


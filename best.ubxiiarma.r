# Função que testa todos os possiveis modelos de uma BARMA e retorna o melhor modelo
# segundo o criterio de AIC.
#
# Implementado por Fabio M Bayer (bayer@ufsm.br) em 02/10/2014
# Adaptado para o modelo UBXII-ARMA por Tatiane F. Ribeiro (tfr1@de.ufpe.br) em 15/10/2020
#
# Dados de entrada:
# - serie: e a serie temporal de interesse;
# - sf:  informacao do inicio e da frequencia da serie temporal que deve estar da seguinte maneira.
#     ex: sf<-c(start=c(1994,7),frequency=12)  # informa o inecio e a frequencia da serie temporal
# - h: quantidade de observacoes que serao reservadas para fazer previsoes dentre deste intervalo da amostra.
# - pmax, qmax, Pmax, Qmax: ordens maximas a serem testadas no modelo SARIMAX. Se nao informar valores serao assumidos os padroes vistos abaixo.

# OBS.: 
#  1- As previsoes para as variaveis explicativas sao feitas com Hol-Winters.
#  2- Essa usa modelos aninhados


best.ubxii<-function(serie, sf, h=6, pmax=6, qmax=6, nbest=10,
                     tau=0.5,link = "logit",X=NA,X_hat=NA)
{
  source("ubxiiarma.fit.r")
  #n<-length(serie)
  #out<-seq((n-h+1),n,1)  # para retirar as ultimas h observacoes
  #y<-ts(serie[-out],start=c(sf[1],sf[2]),frequency=sf[3]) # serie sem as ultimas h observacoes
  y<-ts(serie,start=c(sf[1],sf[2]),frequency=sf[3])
  
  
  # inicializa os criterios AIC
  fit<-ubxiiarma.fit(y, ma=1,diag=0,link = link)
  aicmin<-fit$aic
  
  print(aicmin)
  
  model1<-model2<-model3<-model4<-model5<-0
  model6<-model7<-model8<-model8<-model10<-0
  
  best_aic<-rep(Inf,nbest) # guarda os 10 menores AICs
  #melhores<-rep(0,(nbest)) # guarda as ordem dos 10 melhores modelos
  melhores<-matrix(rep(0,(nbest*3)),ncol=3) # guarda as ordem dos 10 melhores modelos
  colnames(melhores)<-c("p","q","AIC")
  
  tot<-0 # inicializa contador de quantos modelos serao testados
  bug<-0 # incializa contador de quantas vezes deu bug na estimacao Arima(...)
  
  for(p in 0:pmax)
  { 
    for(q in 0:qmax)
    { 
      if(p==0)   ar1<-NA else ar1<-1:p
      if(q==0)   ma1<-NA else ma1<-1:q
      
      if(sum(is.na(c(ar1,ma1)))<2 )
      {   
        print(c(p,q),quote=F)
        fitubxii<-ubxiiarma.fit(y, ar=ar1, ma=ma1,tau=tau,link = link, X=X, X_hat=X_hat)
        tot<-tot+1
        
        if(fitubxii$conv != 0)
        {  
          print(c("NO CONVERGENCE  ",p,q),quote=F)
          bug<-bug+1
          next # sai desse loop e vai para o proximo
        }          
        if(aicmin>fitubxii$aic) # melhor modelo segundo o AIC
        {  
          aicmin<-fitubxii$aic
          best_model_aic <- fitubxii$model
          print("###########################################")
          print(aicmin)
        }
        if(fitubxii$aic<max(best_aic))
        {
          maximo<-order(best_aic)[nbest]
          best_aic[maximo]<-fitubxii$aic
          melhores[maximo,]<-c(p,q,fitubxii$aic)
          #print(melhores)
        }
        
        
      }
    }
  }
  
  
  print(" ",quote=F)
  print("MODELO SELECIONADO VIA AIC",quote=F)
  print(best_model_aic,quote=F)
  print(" ",quote=F)
  print(c("Total de modelos testados =",tot),quote=F)
  print(c("Total erros na estimacao =",bug),quote=F)
  
  print(" ",quote=F)
  print("OS MELHORES MODELOS",quote=F)
  print(melhores,quote=F)
  
}


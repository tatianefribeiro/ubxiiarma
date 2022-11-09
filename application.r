###############################################################################
#Paper: Unit Burr XII autoregressive moving average models
#Tittle: Application of the UBXII-ARMA(p,q) to the actual data
#Author: Tatiane F. Ribeiro
#Last update: October 26, 2022
###############################################################################
rm(list = objects())
# Required packages
library(gdata) #allows to read excel files.  
library(e1071) #allows to calculate asymetry and kustosis.
library(tidyverse)
library(forecast)

# Useful functions and its sources
# ßARMA - available sources at
# https://github.com/vscher/barma           (1 and 2)
# http://www.ufsm.br/bayer/boot-barma.zip   (3)
source("barma.fit.r")    #1
source("barma.r")        #2
source("best.barma.r")   #3

# KARMA: available sources at
# https://github.com/fabiobayer/KARMA
source("kum-mu-phi.r")
source("karma.fit.r")
source("karma.r")
source("best.karma.r")
#source("best.karma_cristiane.r")

# UBXIIARMA
source("ubxiiarma.fit.r")
source("best.ubxiiarma.r")

# Data set
data = read_csv("EAR_southeast_may01-2000_aug31-2019.csv") 

#View(data)

vr <- data$val_eaconsimp4

#enable to accesses variable giving their names.
#attach(data)  

#transform variable vr (stored energy rate) for unit interval.
y1 <- vr/100

month = 5
year = 2000
tau = 0.5

#convert data in time series object.
Y <- ts(y1,
        start = c(year,month),frequency = 12)

#sample size (all observations)
n_first <- length(Y)

#number of forecast steps.
h1 <- 10
month_h = end(Y)[2]-h1+1

#taking off the last h1 observations.
n <- n_first-h1                          #sample size for estimation
y <- ts(y1[1:n],start = c(year,month),frequency = 12)  #in-sample
yh <- ts(Y[(n+1):n_first],start=c(end(Y)[1],month_h), frequency=12) # out-of-sample

hist(y)
monthplot(y)
ggseasonplot(y, year.labels=TRUE, year.labels.left=TRUE) +
  ylab("EAR(%)") +
  ggtitle("Seasonal plot")
#Data description.
summary(y)  #resume
var(y)      #variance.
skewness(y) #skewness.
kurtosis(y) #kurtosis.

#Some graphics, part I
#autocorrelation function.
Acf(y) 

#partial autocorrelation function.
Pacf(y)

#%%%%%%%%%%%%%% WITH REGRESSORS  ########################
# deterministic seazonality
t <- 1:length(y) # in sample
t_hat <- (n+1):(n+h1) # out of sample

C <- cos(2*pi*t/12)
C_hat<-cos(2*pi*t_hat/12) # out of sample

S <- sin(2*pi*t/12)
S_hat<-sin(2*pi*t_hat/12) # out of sample

X <- cbind(C,S)
X_hat <- cbind(C_hat,S_hat)

C_mat<-as.matrix(cos(2*pi*t/12)) # in sample 
S_mat<-as.matrix(sin(2*pi*t/12)) # in sample 

tsbl <- as_tibble(y)

as.integer(time(y)) 

  
which(as.integer(time(y))<2002)  
which(as.integer(time(y))>=2013)
#crisis <- c(rep(1,20),rep(0,132),rep(1,70))
# Lane:dummy para as primeiras 20 obs (8 em 2000 e 12 em 2001) e pras últimas 70=222-152 (obs 153 até 222, tudo pós 2012)
#D <- c(rep(1,20),rep(0,length(y)-20-70),rep(1,70)) #Dummy down
#D_hat <- rep(1, h1)

l1 <- length(which(as.integer(time(y))<2002))
l2 <- which(as.integer(time(y))>=2013)[1]
D <- c(rep(1,l1),
       rep(0,(length(y)-l1-(length(y)-l2+1))),
       rep(1,length(y)-l2+1)) #Dummy down
D_hat <- rep(1, h1)

X_hat_mat <- cbind(as.matrix(C_hat),as.matrix(S_hat), rep(0,h1))
X_mat <- cbind(C_mat,S_mat,D)

X <- cbind(C,S,D)
X_hat <- cbind(C_hat,S_hat,D_hat)

# # Choosing the best models from the UBXII-ARMA, BARMA, and KARMA classes
# pmax = 3
# qmax = 3
# ubxii_best <- best.ubxii(y, sf = c(start = c(year,month),frequency = 12),
#            h=h1, pmax=pmax, qmax=qmax,
#            nbest=10, tau=tau, link = "logit", X = X,X_hat=X_hat)
# 
# barma_best <- best.barma(y, sf = c(start = c(year,month),frequency = 12) ,
#            h=h1, pmax=pmax, qmax=qmax,
#            nbest=10,link = "logit",
#            X = X, X_hat=X_hat)
# 
# karma_best <- best.karma(y, sf = c(start = c(year,month),frequency = 12) ,
#            h=h1, pmax=pmax, qmax=qmax,
#            nbest=10,link = "logit",
#            X = X, X_hat=X_hat)

### SUGESTION FROM REFEREE: TO COMPARE TO THE WHITE NOISE
fit_wn <- arima(y,order=c(2,0,0))
fv_fit_wn <- fitted.values(fit_wn)

plot(y,xlab = "h",ylab = "Stored hydroelectric energy",
     ylim = c(0,1))
lines(fv_fit_wn,lty=2)

summary(fv_fit_wn)
summary(fit_ubxiiarma$fitted)

fit_wn_forec <- arima(yh,c(2,0,0))
prev_fit_wn <- fitted.values(fit_wn_forec)
summary(prev_fit_wn)
summary(fit_ubxiiarma$forecast)
summary(yh)
auto.arima(y)

#The final models are given as follows
p_ubxiiarma <- 1:2#1:2
q_ubxiiarma <- NA
p_barma <- 1:2
q_barma <- NA
p_karma <- 1:2
q_karma <- NA
fit_ubxiiarma <- ubxiiarma.fit(y,ar=p_ubxiiarma,
                               ma=q_ubxiiarma,
                               tau = tau,
                               link = "logit",h=h1, 
                               diag=2, #2:save the plots
                               X = X, X_hat=X_hat)

fit_barma <- barma(y,ar=p_barma,ma=q_barma,link = "logit",
                   h=h1, diag=0, 
                   X = X, X_hat=X_hat)
fit_karma <- karma(y,ar=p_karma,ma=q_karma,link = "logit", 
                   h=h1, diag=0, X = X,
                   X_hat = X_hat)


# Diagnostic measures 
res_UBXII = fit_ubxiiarma$residuals
Box.test(res_UBXII, lag = 10, type =  "Ljung-Box", fitdf = 2)

res_BARMA = fit_barma$resid5  #quant res
Box.test(res_BARMA, lag = 10, type =  "Ljung-Box", fitdf = 2)

res_KARMA = fit_karma$resid3 #quant res
Box.test(res_KARMA, lag = 10, type =  "Ljung-Box", fitdf = 2)

fit_ubxiiarma$aic
fit_barma$aic
fit_karma$aic

fit_ubxiiarma$bic
fit_barma$bic
fit_karma$bic


##################################################################
# Forecasting adequacy measures [max(p,q)]
m_ubxiiarma <- max(p_ubxiiarma,q_ubxiiarma,na.rm = T)
m_barma <- max(p_barma,q_barma,na.rm = T)
m_karma <- max(p_karma,q_karma,na.rm = T)

n1_ubxiiarma <- n-m_ubxiiarma
n1_barma <- n-m_barma
n1_karma <- n-m_karma

##  MSE  ##
mse_ubxiiarma <- ((sum(y[-c(1:m_ubxiiarma)]-fit_ubxiiarma$fitted[(m_ubxiiarma+1):n]))^2)/n1_ubxiiarma
mse_barma <- ((sum(y[-c(1:m_barma)]-fit_barma$fitted[(m_barma+1):n]))^2)/n1_barma
mse_karma <- ((sum(y[-c(1:m_karma)]-fit_karma$fitted[(m_karma+1):n]))^2)/n1_karma

mse_prev_ubxiiarma <- ((sum(yh-fit_ubxiiarma$forecast))^2)/h1
mse_prev_barma <- ((sum(yh-fit_barma$forecast))^2)/h1
mse_prev_karma <- ((sum(yh-fit_karma$forecast))^2)/h1
#WHITE NOISE
mse_prev_wn <- ((sum(yh-prev_fit_wn))^2)/h1

##MAPE##
mape_ubxiiarma <- (sum(abs(y[-c(1:m_ubxiiarma)]-fit_ubxiiarma$fitted[(m_ubxiiarma+1):n])/abs(y[-c(1:m_ubxiiarma)])))/n1_ubxiiarma
mape_barma <- (sum(abs(y[-c(1:m_barma)]-fit_barma$fitted[(m_barma+1):n])/abs(y[-c(1:m_barma)])))/n1_barma
mape_karma <- (sum(abs(y[-c(1:m_karma)]-fit_karma$fitted[(m_karma+1):n])/abs(y[-c(1:m_karma)])))/n1_ubxiiarma

mape_prev_ubxiiarma<-(sum(abs(yh-fit_ubxiiarma$forecast)/abs(yh)))/h1 
mape_prev_barma<-(sum(abs(yh-fit_barma$forecast)/abs(yh)))/h1 
mape_prev_karma<-(sum(abs(yh-fit_karma$forecast)/abs(yh)))/h1 
#WHITE NOISE
mape_prev_wn <-(sum(abs(yh-prev_fit_wn)/abs(yh)))/h1 

##MASE##
mase_ubxiiarma <- (1/n1_ubxiiarma)*(sum(abs(y[-c(1:m_ubxiiarma)]-fit_ubxiiarma$fitted[(m_ubxiiarma+1):n])/
                                          ((1/(n1_ubxiiarma-1))*sum(abs(diff(y[-c(1:m_ubxiiarma)]))))))
mase_barma <- (1/n1_barma)*(sum(abs(y[-c(1:m_barma)]-fit_barma$fitted[(m_barma+1):n])/
                                  ((1/(n1_barma-1))*sum(abs(diff(y[-c(1:m_barma)]))))))
mase_karma <- (1/n1_karma)*(sum(abs(y[-c(1:m_karma)]-fit_karma$fitted[(m_karma+1):n])/
                                  ((1/(n1_karma-1))*sum(abs(diff(y[-c(1:m_karma)]))))))

mase_prev_ubxiiarma <- (1/h1)*(sum(abs(yh-fit_ubxiiarma$forecast)/
                                     ((1/(h1-1))*sum(abs(diff(yh))))))
mase_prev_barma <- (1/h1)*(sum(abs(yh-fit_barma$forecast)/
                                 ((1/(h1-1))*sum(abs(diff(yh))))))
mase_prev_karma <- (1/h1)*(sum(abs(yh-fit_karma$forecast)/
                                 ((1/(h1-1))*sum(abs(diff(yh))))))
#WHITE NOISE
mase_prev_wn <- (1/h1)*(sum(abs(yh-prev_fit_wn)/
                                 ((1/(h1-1))*sum(abs(diff(yh))))))


###Results###
res_forec_meas <- matrix(NA, nrow = 6,ncol = 3)
rownames(res_forec_meas) <- c("MSE (in)", "MAPE (in)", "MASE (in)",
                              "MSE (out)", "MAPE (out)", "MASE (out)")
colnames(res_forec_meas) <- c("UBXII-ARMA", "BARMA", "KARMA")
res_forec_meas[,1] <- c(mse_ubxiiarma,
                        mape_ubxiiarma,
                        mase_ubxiiarma,
                        mse_prev_ubxiiarma,
                        mape_prev_ubxiiarma,
                        mase_prev_ubxiiarma
)

res_forec_meas[,2] <- c(mse_barma,
                        mape_barma,
                        mase_barma,
                        mse_prev_barma,
                        mape_prev_barma,
                        mase_prev_barma
)
res_forec_meas[,3] <- c(mse_karma,
                        mape_karma,
                        mase_karma,
                        mse_prev_karma,
                        mape_prev_karma,
                        mase_prev_karma
)
round(res_forec_meas,4)

stargazer::stargazer(res_forec_meas,digits=4)

#Final Fitted models
rbind(fit_ubxiiarma$model,fit_barma$model,fit_karma$model)
stargazer::stargazer(rbind(fit_ubxiiarma$model,fit_barma$model,fit_karma$model), digits = 4)

#stargazer::stargazer(summary(y),digits = 4)  #resume
var(y)      #variance.
skewness(y) #skewness.
kurtosis(y) #kurtosis.


w1<-5
h2<-4
postscript(file = "forecast_comp.eps",horizontal=F,paper="special",width = w1, height = h2,family = "Times")
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  Prev_val<-as.numeric(yh)
  pred_val_UBXII<-fit_ubxiiarma$forecast
  pred_val_BETA<-fit_barma$forecast
  pred_val_KARMA<-fit_karma$forecast
  plot(Prev_val,xlab = "h",ylab = "Stored energy proportion",
       ylim = c(0,1))
  lines(pred_val_UBXII,lty=1)
  lines(pred_val_BETA,lty=2)
  lines(pred_val_KARMA,lty=3)
  #lines(as.numeric(prev_fit_wn), lty=4,col=2)
  legend("topleft",c("Observed data","UBXII-AR(2)",
                     expression(paste(beta,"AR(2)"), "KAR(2)") ),
         #pch=vpch,
         pt.bg="white", pch = c(1,NA,NA,NA),lty = c(NA,1,2,3), bty="n")
  
}
dev.off()


w1<-5
h2<-4
postscript(file = "HW.eps",horizontal=F,paper="special",
           width = w1, height = h2,family = "Times")
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
AE <- HoltWinters(y,seasonal = "multiplicative",)
AE
pred_val_AE<- as.numeric(predict(AE,h1))
Prev_val<-as.numeric(yh)
pred_val_UBXII<-fit_ubxiiarma$forecast
pred_val_BETA<-fit_barma$forecast
pred_val_KARMA<-fit_karma$forecast
plot(Prev_val,xlab = "h",ylab = "Stored hydroelectric energy",
     ylim = c(0,1))
lines(pred_val_UBXII,lty=1)
lines(pred_val_BETA,lty=2)
lines(pred_val_KARMA,lty=3)
lines(pred_val_AE,lty=4,col=2)
legend("bottomleft",c("Actual values","UBXII-AR(2)",
                   expression(paste(beta,"AR(2)"), "KAR(2)", 
                              "HOLT-Winters" ) ),#pch=vpch,
       pt.bg="white", pch = c(1,NA,NA,NA,NA),lty = c(NA,1,2,3,4),col = c(1,1,1,1,2), bty="n")
}
dev.off()


w1<-6
h2<-5
postscript(file = "stl.eps",horizontal=F,
           paper="special",width = w1, height = h2,
           family = "Times",
           pointsize = 15)
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
plot(stl(y, s.window=12))
}
dev.off()


w1<-6
h2<-5
postscript(file = "seasonality.eps",horizontal=F,
           paper="special",width = w1, height = h2,
           family = "Times",
           pointsize = 15)
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  monthplot(y,ylab = "Stored hydroelectric energy")
}
dev.off()


w1<-6
h2<-5
postscript(file = "y_acf.eps",horizontal=F,
           paper="special",width = w1, height = h2,
           family = "Times",
           pointsize = 15)
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  Acf(y,lag=30, main="")
}
dev.off()

w1<-6
h2<-5
postscript(file = "y_pacf.eps",horizontal=F,
           paper="special",width = w1, height = h2,
           family = "Times",
           pointsize = 15)
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  Pacf(y, main="")
}
dev.off()




##### Without effect AR terms
alpha_hat <- fit_ubxiiarma$alpha
xbeta <- X%*%fit_ubxiiarma$beta
eta_hat <- alpha_hat+xbeta

q_hat <- ts(as.vector(exp(eta_hat)/(1+exp(eta_hat))))

w1<-7
h2<-5
postscript(file = "fitt_value.eps",horizontal=F,
           paper="special",width = w1, height = h2,
           family = "Times",
           pointsize = 15)
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  plot(y,type="l",ylab="Stored hydroelectric energy",xlab="Time",
       ylim=c(min(y),max(y)))
  lines(fit_ubxiiarma$fitted,col="blue",lty=2)
  lines(ts(q_hat, start=c(2000,7), freq=12), col="red", 
        lty=2) #Lane
  legend("topright",c("Observed data","Predicted median","Naive prediction"),#pch=vpch,
         pt.bg="white", lty=c(1,2,2), bty="n",col=c(1,"blue",
                                                    "red"))
  
}
dev.off()


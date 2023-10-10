##############################################################
# Paper: Forecasting the proportion of stored energy
#        using the Unit Burr XII quantile
#        autoregressive moving average model (UBXII-ARMA)
# Author: Tatiane F. Ribeiro
# Last update: October 06, 2023
#############################################################
# Clear the memory
rm(list = objects())

# Required packages 
library(e1071)     # To calculate asymmetry and kurtosis
library(tidyverse) # To use the read_csv() function

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

# UBXIIARMA
source("ubxiiarma.fit.r")
source("best.ubxiiarma.r")

# Data set
data = read_csv("EAR_southeast_may01-2000_aug31-2019.csv") 
# View(data)

# Proportion of stored energy 
y1 = data$val_eaconsimp4/100

month = 5
year  = 2000
tau   = 0.5

# Convert data in time series object
Y = ts(y1,
       start     = c(year,month),
       frequency = 12)

# Sample size (all observations)
n_first = length(Y)

# Number of forecast steps
h1      = 10
month_h = end(Y)[2]-h1+1

# Taking off the last h1 observations
n = n_first-h1  
y = ts(y1[1:n],          # Training data set
       start     = c(year,month),
       frequency = 12) 

yh = ts(Y[(n+1):n_first], # Test data set
        start     = c(end(Y)[1],month_h), 
        frequency = 12) 
summary(yh)

# Figure 2: Some plots of y
plot(stl(y, s.window = 12))
plot(y)

# Seasonality plots
monthplot(y)

# Autocorrelation function.
acf(y) 

# Partial autocorrelation function.
pacf(y)

# Table 6: Descriptive statistics of y 
summary(y)  # resume
var(y)      # variance
skewness(y) # skewness
kurtosis(y) # kurtosis

# REGRESSORS 
# Deterministic seasonality
t     = 1:length(y)      # in-sample
t_hat = (n+1):(n+h1)     # out-of-sample

C     = cos(2*pi*t/12)
C_hat = cos(2*pi*t_hat/12)  

S     = sin(2*pi*t/12)
S_hat = sin(2*pi*t_hat/12) 

# Defining the dummy variable for the water crisis
# crisis = c(rep(1,20),rep(0,132),rep(1,70))
# dummy = 1, for  t = 1,...,20:
# 8 obs in 2000 and 12 in 2001
# and dummy = 1, for last 70 obs:
# 70 = 222-152 (obs 153 until 222, all after 2012)

which(as.integer(time(y)) < 2002)  
which(as.integer(time(y)) >= 2013)

l1    = length(which(as.integer(time(y)) < 2002))
l2    = which(as.integer(time(y)) >= 2013)[1]
D     = c(rep(1, l1),
          rep(0, (length(y)-l1-(length(y)-l2+1))),
          rep(1, length(y)-l2+1)) 
D_hat = rep(1, h1)

# Regressors matrix
X     = cbind(C,S,D)                # in-sample
X_hat = cbind(C_hat,S_hat,D_hat)    # out-of-sample

# Choosing the best models from the UBXII-ARMA, BARMA, and KARMA classes
pmax = 3
qmax = 3
ubxii_best = best.ubxii(y, sf = c(start = c(year,month),frequency = 12),
           h = h1, pmax = pmax, qmax = qmax,
           nbest = 10, tau = tau, link = "logit", X = X, X_hat = X_hat)

barma_best = best.barma(y, sf = c(start = c(year,month), frequency = 12) ,
           h = h1, pmax = pmax, qmax = qmax,
           nbest = 10,link = "logit",
           X = X, X_hat = X_hat)

karma_best = best.karma(y, sf = c(start = c(year,month),frequency = 12) ,
           h = h1, pmax = pmax, qmax = qmax,
           nbest = 10,link = "logit",
           X = X, X_hat = X_hat)

# The final models are given as follows
p_ubxiiarma = 1:2
q_ubxiiarma = NA
p_barma = 1:2
q_barma = NA
p_karma = 1:2
q_karma = NA
fit_ubxiiarma = ubxiiarma.fit(y,
                              ar    = p_ubxiiarma,
                              ma    = q_ubxiiarma,
                              tau   = tau,
                              link  = "logit",
                              h     = h1, 
                              diag  = 1,
                              X     = X, 
                              X_hat = X_hat)

fit_barma = barma(y,
                  ar    = p_barma,
                  ma    = q_barma,
                  link  = "logit",
                  h     = h1, 
                  diag  = 1,
                  X     = X,
                  X_hat = X_hat)

fit_karma = karma(y,
                  ar    = p_karma,
                  ma    = q_karma,
                  link  = "logit", 
                  h     = h1,
                  diag  = 1,
                  X     = X,
                  X_hat = X_hat)

# Table 7: Final fitted models 
rbind(fit_ubxiiarma$model,fit_barma$model,fit_karma$model)
#stargazer::stargazer(rbind(fit_ubxiiarma$model,fit_barma$model,fit_karma$model), digits = 4)

# Residual ACF and PACF
acf(fit_ubxiiarma$residuals)
pacf(fit_ubxiiarma$residuals)
acf(fit_barma$resid5)
pacf(fit_barma$resid5)
acf(fit_karma$resid3)
pacf(fit_karma$resid3)

# Out-of-sample forecasting (computation of the quantities from Table 8)
# Auxiliary matrix
X_hat1 = cbind(C_hat,S_hat,D_hat)
# Test sets
test_data  = cbind(yh,X_hat)

# Create a vector to store the forecasted values
forecasts_ubxiiarma = forecasts_karma = forecasts_barma = matrix(NA, nrow = 10, ncol = 10)
diag = 0
for (i in 1:nrow(test_data)) {
  
  # Current training set
  if(i == 1){
    X_hat = t(as.matrix(X_hat1[1:i,]))
  }else{
    X_hat = X_hat1[1:i,]
  }
  
  # Fit UBXII-AR(2) model to the current training set
  fit_ubxiiarma = ubxiiarma.fit(y,
                                ar    = 1:2,
                                ma    = NA,
                                tau   = tau,
                                link  = "logit",
                                h     = i, 
                                diag  = diag, 
                                X     = X,
                                X_hat = X_hat)
  
  # Fit betaAR(2) model to the current training set
  fit_barma = barma(y,
                    ar    = 1:2,
                    ma    = NA,
                    link  = "logit",
                    h     = i, 
                    diag  = diag, 
                    X     = X,
                    X_hat = X_hat)
  
  # Fit KAR(2) model to the current training set
  fit_karma = karma(y,
                    ar    = 1:2,
                    ma    = NA,
                    link  = "logit",
                    h     = i, 
                    diag  = diag, 
                    X     = X,
                    X_hat = X_hat)
  
  # Store the forecasted values
  forecasts_ubxiiarma[i,1:i] = fit_ubxiiarma$forecast
  forecasts_barma[i,1:i]     = fit_barma$forecast
  forecasts_karma[i,1:i]     = fit_karma$forecast
}

# Print the forecasted values
print("Forecasted values (UBXII-AR(2)):")
round(forecasts_ubxiiarma,4)
print("Forecasted values (betaAR(2)):")
round(forecasts_barma,4)
print("Forecasted values (KAR(2)):")
round(forecasts_karma,4)

# Null vectors to salve the MSEs and MAPEs 
mse_ubxiiarma = mape_ubxiiarma = numeric(h1)
mse_barma     = mape_barma     = numeric(h1)
mse_karma     = mape_karma     = numeric(h1)
for (i in 1:h1) {
  # Calculate the actual values from the test set
  actual_values = yh[1:i]
  
  forecasts_current_ubxiiarma = as.numeric(na.omit(forecasts_ubxiiarma[i,]))
  forecasts_current_barma     = as.numeric(na.omit(forecasts_barma[i,]))
  forecasts_current_karma     = as.numeric(na.omit(forecasts_karma[i,]))
  
  # Calculate evaluation metrics
  mse_ubxiiarma[i]  = mean((forecasts_current_ubxiiarma - actual_values)^2)
  mape_ubxiiarma[i] = mean(abs((forecasts_current_ubxiiarma - actual_values) / actual_values)) * 100
  
  mse_barma[i]  = mean((forecasts_current_barma - actual_values)^2)
  mape_barma[i] = mean(abs((forecasts_current_barma - actual_values) / actual_values)) * 100
  
  mse_karma[i]  = mean((forecasts_current_karma - actual_values)^2)
  mape_karma[i] = mean(abs((forecasts_current_karma - actual_values) / actual_values)) * 100
}

# Evaluation metrics
results = rbind(mse_ubxiiarma,
                mape_ubxiiarma,
                mse_barma,
                mape_barma,
                mse_karma,
                mape_karma
)

# Table 8
round(results, 4)
# stargazer::stargazer(results, digits = 4)

# Saving the plots (Figure 4a)
# Without effect AR terms
alpha_hat = fit_ubxiiarma$alpha
xbeta = X%*%fit_ubxiiarma$beta
eta_hat = alpha_hat+xbeta
q_hat = ts(as.vector(exp(eta_hat)/(1+exp(eta_hat))))

# Plot dimension
w1 = 7
h2 = 5

# Figure 4a
# Save a eps file
postscript(file       = "fitt_value.eps",
           horizontal = F,
           paper      = "special",
           width      = w1, 
           height     = h2,
           family     = "Times",
           pointsize  = 15)
{
  par(mfrow = c(1,1))
  par(mar = c(2.8, 2.7, 1, 1))  
  par(mgp = c(1.7, 0.45, 0))
  plot(y,
       type = "l",
       ylab = "Stored hydroelectric energy",
       xlab = "Time",
       ylim = c(min(y),max(y)))
  lines(fit_ubxiiarma$fitted,
        col = "blue",
        lty = 2)
  lines(ts(q_hat, 
           start = c(year, month), 
           freq = 12), 
        col = "red", 
        lty = 2)
  legend("topright",
         c("Observed data","Predicted median","Naive prediction"),
         pt.bg = "white", 
         lty = c(1,2,2), 
         bty = "n",
         col = c(1,"blue","red"))
}
dev.off()




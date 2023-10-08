###############################################################################
# Paper: Forecasting the proportion of stored energy
#        using the Unit Burr XII quantile
#        autoregressive moving average model (UBXII-ARMA)
# Author: Tatiane F. Ribeiro
# Last update: October 07, 2023
###############################################################################
# Clear the memory
rm(list = objects())

# Required packages
library(tidyverse) # To use the read_csv() function

# UBXII-ARMA with r_t: quantile residuals
source("ubxiiarma.fit_qr.r")
source("best.ubxiiarma_qr.r")

# Data set
data = read_csv("EAR_southeast_may01-2000_aug31-2019.csv") 
# View(data)

# Proportion of stored energy 
y1 = data$val_eaconsimp4/100

month = 5
year  = 2000
tau   = 0.9

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

#############################
# Deterministic seazonality
#############################
t     = 1:length(y) # in-sample
t_hat = (n+1):(n+h1) # out-of-sample

C     = cos(2*pi*t/12)
C_hat = cos(2*pi*t_hat/12)  

S     = sin(2*pi*t/12)
S_hat = sin(2*pi*t_hat/12)  

j      = 2
C2     = cos(2*pi*t*j/12)
C2_hat = cos(2*pi*t_hat*j/12)  

S2     = sin(2*pi*t*j/12)
S2_hat = sin(2*pi*t_hat*j/12)  

j      = 4
C4     = cos(2*pi*t*j/12)
C4_hat = cos(2*pi*t_hat*j/12)  

S4     = sin(2*pi*t*j/12)
S4_hat = sin(2*pi*t_hat*j/12)  

X     = cbind(C,C2, S4)
X_hat = cbind(C_hat, C2_hat, S4_hat)

# Choosing the best UBXII-ARMA model according to AIC criterion
pmax = 3
qmax = 3
ubxii_best = best.ubxii(y, sf = c(start = c(year,month),frequency = 12),
           h = h1, pmax = pmax, qmax = qmax,
           nbest = 10, tau = tau, link = "logit", X = X, X_hat = X_hat)

# Fit of the selected model
p_ubxiiarma = 1:2
q_ubxiiarma = 1:3

fit_ubxiiarma = ubxiiarma.fit(y,
                              ar    = p_ubxiiarma,
                              ma    = q_ubxiiarma,
                              tau   = tau,
                              link  = "logit",
                              h     = h1, 
                              X     = X, 
                              X_hat = X_hat)

fit_ubxiiarma$model

# Residuals
res_UBXII = fit_ubxiiarma$residuals
acf(res_UBXII)
pacf(res_UBXII)

# Plots for the paper
# Fitted against observed values
w1 = 7
h2 = 5
postscript(file       = "fitt_value_tau090.eps",
           horizontal = F,
           paper      = "special",
           width      = w1, 
           height     = h2,
           family     = "Times",
           pointsize  = 15)
{
  par(mfrow = c(1,1))
  par(mar   = c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp   = c(1.7, 0.45, 0))
  
  plot(y,
       type = "l",
       ylab = "Stored hydroelectric energy",
       xlab = "Time",
       ylim = c(min(y),max(y)))
  
  lines(fit_ubxiiarma$fitted,
        col = "blue",
        lty = 2)
  legend("topright",c("Observed data",
                      "Pred. 90th perc."),
         pt.bg = "white",
         lty = c(1,2), 
         bty = "n",
         col = c(1,"blue"))
  
}
dev.off()

# UBXII-ARMA residuals acf
w1 = 6
h2 = 5
postscript(file       = "res_acf_tau090.eps",
           horizontal = F,
           paper      = "special",
           width      = w1, 
           height     = h2,
           family     = "Times",
           pointsize  = 15)
{
  par(mfrow = c(1,1))
  par(mar   = c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp   = c(1.7, 0.45, 0))
  acf(res_UBXII, 
      lag  = 30, 
      main = "")
}
dev.off()

# UBXII-ARMA residuals pacf
postscript(file       = "res_pacf_tau090.eps",
           horizontal = F,
           paper      = "special",
           width      = w1, 
           height     = h2,
           family     = "Times",
           pointsize  = 15)
{
  par(mfrow = c(1,1))
  par(mar   = c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp   = c(1.7, 0.45, 0))
  pacf(res_UBXII, 
      lag  = 30, 
      main = "")
}
dev.off()




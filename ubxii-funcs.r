# cumulative distribution function
cdf_UBXII <- function(y,q,c,tau=0.5)   #UBXII cdf (useful for qqplot)
{
  (1+(log(1/y))^c)^(log(tau)/log(1+(log(1/q))^c))
}


# inversion method for randon generation
r_UBXII <- function(n,q_t,c,tau=0.5)   #It generates occurences of Y_i ~ UBXII (q_t, c)
{
  u = runif(n)
  y_UBXII = exp(-(u^(-1/(log(1/tau)/log(1+(log(1/q_t))^c)))-1)^(1/c)) #UBXII qf
  return(y_UBXII)
}


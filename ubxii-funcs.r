# cumulative distribution function
cdf_UBXII <- function(y,q,c,tau=0.5)    
{
  (1+(log(1/y))^c)^(log(tau)/log(1+(log(1/q))^c))
}


# inversion method to generate samples from UBXII distribution (Y_t ~ UBXII (q_t, c))
r_UBXII <- function(n,q_t,c,tau=0.5)   
{
  u = runif(n)
  y_UBXII = exp(-(u^(-1/(log(1/tau)/log(1+(log(1/q_t))^c)))-1)^(1/c)) # UBXII qf
  return(y_UBXII)
}


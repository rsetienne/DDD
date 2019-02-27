dd_loglik_high_lambda <- function(pars1,pars2,brts)
{
  lbrts = length(brts)
  if(brts[lbrts] == 0)
  {
    brts = brts[-lbrts]
    lbrts = length(brts)
  }
  soc = pars2[6]
  N = lbrts + (soc - 1)
  startk = soc
  endk = N
  mu = pars1[2]
  brtsdiff = c(brts[2:lbrts],0) - c(brts[1:lbrts])   
  
  loglik <- log(N) +
    (soc == 2) * (log(N + 1) - log(6)) + 
    (N - 1) * log(mu) +
    (1 + (pars2[4] == 0)) * lgamma(N) +
    - (N - 1) * log(N - 1) +
    - (soc == 2) * log(mu) + 
    - mu/(N - 1) * sum((startk:endk)*((startk - 1):(endk - 1)) * brtsdiff[1:lbrts])
  return(loglik)
}
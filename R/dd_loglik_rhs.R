dd_loglik_rhs_precomp = function(pars,x)
{  
  lx = length(x)
  la = pars[1]
  mu = pars[2]
  K = pars[3]
  if (length(pars) < 6) {
    kk = pars[4]
    ddep = pars[5]
  } else {
    r = pars[4]
    kk = pars[5]
    ddep = pars[6]
    phi = ifelse(r == Inf, 1, r / (1 + r)) # else r/(1+r) can be NaN
    eq_rate = phi * la + (1 - phi) * mu 
  }
  n0 = (ddep == 2 | ddep == 4)
  
  nn <- c(0, 0:(lx + 2 * kk))
  lnn <- length(nn)
  
  if (ddep == 1) {
    lavec = pmax(0, la - (la - mu) / K * nn)
    muvec = rep(mu, lnn)
  } else if (ddep == 1.3) {
    lavec = pmax(0, la * (1 - nn / K))
    muvec = rep(mu, lnn)
  } else if (ddep == 1.4) {
    lavec = pmax(0, la * nn / (nn + K))
  } else if(ddep == 1.5) {
    lavec = pmax(0, la * nn / K * (1 - nn / K))
    muvec = rep(mu, lnn)
  } else if (ddep == 2 | ddep == 2.1 | ddep == 2.2) {
    y = -(log(la / mu) / log(K + n0)) ^ (ddep != 2.2)
    lavec = pmax(0, la * (nn + n0) ^ y)
    muvec = rep(mu, lnn)
  } else if(ddep == 2.3) {
    y = -K
    lavec = pmax(0, la * (nn + n0) ^ y)
    muvec = rep(mu, lnn)
  } else if (ddep == 3) {
    lavec = rep(la, lnn)
    muvec = mu + (la - mu) * nn / K
  } else if (ddep == 4 | ddep == 4.1 | ddep == 4.2) {
    lavec = rep(la, lnn)
    y = (log(la / mu) / log(K + n0)) ^ (ddep != 4.2)
    muvec = mu * (nn + n0) ^ y
  } else if (ddep == 5) {
    lavec = pmax(0, la - (la - eq_rate) * nn / K)
    muvec = mu + (eq_rate - mu) / K * nn
  } else if (ddep == 6) {
    y = log(eq_rate / mu) / log(K)
    lavec = pmax(0, la - (la - eq_rate) * nn / K)
    muvec = mu * nn ^ y
  } else if (ddep == 7) {
    y1 = -log(la / eq_rate) / log(K)
    y2 = log(eq_rate / mu) / log(K)
    lavec = pmax(0, la * nn ^ y1)
    muvec = mu * nn ^ y2
  } else if (ddep == 8) {
    y = -log(la / eq_rate) / log(K)
    lavec = pmax(0, la * nn ^ y)
    muvec = mu + (eq_rate - mu) / K * nn
  } else if (ddep == 9) {
    y = log(la / mu) / K
    lavec = pmax(0, la * exp(-nn * y))
    muvec = rep(mu, lnn)
  } else if (ddep == 10) {
    y = log(la / mu) / K
    lavec = rep(la, lnn)
    muvec = mu * exp(nn * y)
  } else if (ddep == 11) {
    y = log(eq_rate / mu) / K
    lavec = pmax(0, la - (la - eq_rate) * nn / K )
    muvec = mu * exp(nn * y)
  } else if (ddep == 12) {
    y1 = log(la / eq_rate) / K
    y2 = log(eq_rate / mu) / K
    lavec = pmax(0, la * exp(-nn * y1))
    muvec = mu * exp(nn * y2)
  } else if (ddep == 13) {
    y = log(la / eq_rate) / K
    lavec = pmax(0, la * exp(-nn * y))
    muvec = mu + (eq_rate - mu) / K * nn
  } else if (ddep == 14) {
    y1 = log(la / eq_rate) / K
    y2 = log(eq_rate / mu) / log(K)
    lavec = pmax(0, la * exp(-nn * y1))
    muvec = mu * nn ^ y2
  } else if (ddep == 15) {
    y1 = -log(la / eq_rate) / log(K)
    y2 = log(eq_rate / mu) / K
    lavec = pmax(0, la * nn ^ y1)
    muvec = mu * exp(nn * y2)
  }
  return(c(lavec, muvec, nn))
}  

dd_loglik_rhs = function(t,x,parsvec)
{
  # Unpack parameter vector
  lv = (length(parsvec) - 1)/3
  lavec = parsvec[1:lv]
  muvec = parsvec[(lv + 1):(2 * lv)]
  nn = parsvec[(2 * lv + 1):(3 * lv)]
  k = parsvec[length(parsvec)] # nb species in phylo at time t
  
  lx = length(x)
  qn_vec = c(0, x, 0)
  nvec <- 2:(lx + 1)
  # DDD master system
  dx <-  lavec[nvec + k - 1] * nn[nvec + 2 * k - 1] * qn_vec[nvec - 1] + 
    muvec[nvec + k + 1] * nn[nvec + 1] * qn_vec[nvec + 1] - 
    (lavec[nvec + k] + muvec[nvec + k]) * nn[nvec + k] * qn_vec[nvec]
  return(list(dx))
}
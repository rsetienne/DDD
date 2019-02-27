dd_loglik_rhs_precomp = function(pars,x)
{  
  lx = length(x)
  la = pars[1]
  mu = pars[2]
  K = pars[3]
  if(length(pars) < 6)
  {
    kk = pars[4]
    ddep = pars[5]
  } else {
    r = pars[4]
    kk = pars[5]
    ddep = pars[6]
  }
  n0 = (ddep == 2 | ddep == 4)
  
  nn = -1:(lx+2*kk)
  lnn = length(nn)
  nn = pmax(rep(0,lnn),nn)
  
  if(ddep == 1)
  {
    lavec = pmax(rep(0,lnn),la - (la - mu)/K * nn)
    muvec = mu * rep(1,lnn)
  } else {
    if(ddep == 1.3)
    {
      lavec = pmax(rep(0,lnn),la * (1 - nn/K))
      muvec = mu * rep(1,lnn)
    } else {
      if(ddep == 1.4)
      {
        lavec = pmax(rep(0,lnn),la * nn/(nn + K))
      } else {
        if(ddep == 1.5)
        {
          lavec = pmax(rep(0,lnn),la * nn/K * (1 - nn/K))
          muvec = mu * rep(1,lnn)
        } else {
          if(ddep == 2 | ddep == 2.1 | ddep == 2.2)
          {
            y = -(log(la/mu)/log(K+n0))^(ddep != 2.2)
            lavec = pmax(rep(0,lnn),la * (nn + n0)^y)
            muvec = mu * rep(1,lnn)
          } else {
            if(ddep == 2.3)
            {
              y = K
              lavec = pmax(rep(0,lnn),la * (nn + n0)^y)
              muvec = mu * rep(1,lnn)
            } else {
              if(ddep == 3)
              {
                lavec = la * rep(1,lnn)
                muvec = mu + (la - mu) * nn/K
              } else {
                if(ddep == 4 | ddep == 4.1 | ddep == 4.2)
                {
                  lavec = la * rep(1,lnn)
                  y = (log(la/mu)/log(K+n0))^(ddep != 4.2)
                  muvec = mu * (nn + n0)^y
                } else {
                  if(ddep == 5)
                  { 
                    lavec = pmax(rep(0,lnn),la - 1/(r + 1)*(la - mu)/K * nn)
                    muvec = muvec = mu + r/(r + 1)*(la - mu)/K * nn
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(c(lavec,muvec,nn))
}  

dd_loglik_rhs = function(t,x,parsvec)
{
  lv = (length(parsvec) - 1)/3
  lavec = parsvec[1:lv]
  muvec = parsvec[(lv + 1):(2 * lv)]
  nn = parsvec[(2 * lv + 1):(3 * lv)]
  kk = parsvec[length(parsvec)]
  lx = length(x)
  xx = c(0,x,0)
  dx = lavec[(2:(lx+1))+kk-1] * nn[(2:(lx+1))+2*kk-1] * xx[(2:(lx+1))-1] + muvec[(2:(lx+1))+kk+1] * nn[(2:(lx+1))+1] * xx[(2:(lx+1))+1] - (lavec[(2:(lx+1))+kk] + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]
  return(list(dx))
}
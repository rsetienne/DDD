dd_loglik_bw_rhs = function(t,x,parsvec)
{
  #parsvec = c(dd_loglik_rhs_precomp(pars,x),pars[length(pars) - 1])
  lv = (length(parsvec) - 1)/3
  lavec = parsvec[1:lv]
  muvec = parsvec[(lv + 1):(2 * lv)]
  nn = parsvec[(2 * lv + 1):(3 * lv)]
  kk = parsvec[length(parsvec)]
  lx = length(x) - 1

  xx = c(0,x[1:lx],0)
  
  dx = lavec[(2:(lx+1))+kk] * nn[(2:(lx+1))+2*kk] * xx[(2:(lx+1))+1] + muvec[(2:(lx+1))+kk] * nn[(2:(lx+1))] * xx[(2:(lx+1))-1] - (c(lavec[(2:(lx))+kk],0) + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]
  #dx = lavec[(2:(lx+1))+kk] * nn[(2:(lx+1))+2*kk] * xx[(2:(lx+1))+1] + muvec[(2:(lx+1))+kk] * nn[(2:(lx+1))] * xx[(2:(lx+1))-1] - (lavec[(2:(lx+1))+kk] + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]
  
  dG = x[1 + (kk == 0)]
  
  return(list(c(dx,dG)))
}
# brts = branching times (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# - pars1[1] = la = (initial) speciation rate
# - pars1[2] = mu = extinction rate
# - pars1[3] = K = carrying capacity
# - pars1[4] = r = ratio of diversity-dependence in extinction rate over speciation rate
# - pars2[1] = lx = length of ODE variable x
# - pars2[2] = ddep = diversity-dependent model,mode of diversity-dependence
#  . ddep == 1 : linear dependence in speciation rate with parameter K
#  . ddep == 1.3 : linear dependence in speciation rate with parameter K'
#  . ddep == 1.4 : positive and negative linear diversity-dependence in speciation rate with parameter K'
#  . ddep == 2 : exponential dependence in speciation rate
#  . ddep == 2.1: variant with offset at infinity
#  . ddep == 2.2: 1/n dependence in speciation rate
#  . ddep == 2.3: exponential dependence in speciation rate with parameter x
#  . ddep == 3 : linear dependence in extinction rate
#  . ddep == 4 : exponential dependence in extinction rate
#  . ddep == 4.1: variant with offset at infinity
#  . ddep == 4.2: 1/n dependence in speciation rate
#  . ddep == 5 : linear dependence in speciation and extinction rate
# - pars2[3] = cond = conditioning
#  . cond == 0 : conditioning on stem or crown age
#  . cond == 1 : conditioning on stem or crown age and non-extinction of the phylogeny
#  . cond == 2 : conditioning on stem or crown age and on the total number of extant taxa (including missing species)
#  . cond == 3 : conditioning on the total number of extant taxa (including missing species)
# - pars2[4] = btorph = likelihood of branching times (0) or phylogeny (1), differ by a factor (S - 1)! where S is the number of extant species
# - pars2[5] = parameters and likelihood should be printed (1) or not (0)
# - pars2[6] = likelihood is for a tree with crown age (2) or stem age (1)
# missnumspec = number of missing species    
# methode = the method used in the numerical solving of the set of the ode's or 'analytical' which means matrix exponentiation is used

dd_loglik_test = function(pars1,pars2,brts,missnumspec,methode = 'analytical',rhs_func_name = 'dd_loglik_rhs')
{
  if(methode == 'analytical')
  {
    out = dd_loglik2(pars1,pars2,brts,missnumspec)
  } else
  {
    out = dd_loglik1(pars1,pars2,brts,missnumspec,methode = methode,rhs_func_name = rhs_func_name)
  }
  return(out)
}



#' Loglikelihood for diversity-dependent diversification models
#' 
#' This function computes loglikelihood of a diversity-dependent
#' diversification model for a given set of branching times and parameter
#' values.
#' 
#' 
#' @param pars1 Vector of parameters: 
#' \cr \cr \code{pars1[1]} corresponds to
#' lambda (speciation rate)
#' \cr \code{pars1[2]} corresponds to mu (extinction
#' rate)
#' \cr \code{pars1[3]} corresponds to K (clade-level carrying capacity)
#' @param pars2 Vector of model settings: 
#' \cr \cr \code{pars2[1]} sets the
#' maximum number of species for which a probability must be computed.  This
#' must be larger than 1 + missnumspec + length(brts).
#' \cr \cr \code{pars2[2]}
#' sets the model of diversity-dependence:
#' \cr - \code{pars2[2] == 1} linear
#' dependence in speciation rate with parameter K (= diversity where speciation
#' = extinction)
#' \cr - \code{pars2[2] == 1.3} linear dependence in speciation
#' rate with parameter K' (= diversity where speciation = 0)
#' \cr - \code{pars2[2] == 1.4} : positive diversity-dependence in speciation rate
#' with parameter K' (= diversity where speciation rate reaches half its
#' maximum); lambda = lambda0 * S/(S + K') where S is species richness
#' \cr - \code{pars2[2] == 1.5} : positive and negative diversity-dependence in
#' speciation rate with parameter K' (= diversity where speciation = 0); lambda
#' = lambda0 * S/K' * (1 - S/K') where S is species richness
#' \cr - \code{pars2[2] == 2} exponential dependence in speciation rate with
#' parameter K (= diversity where speciation = extinction)
#' \cr - \code{pars2[2]
#' == 2.1} variant of exponential dependence in speciation rate with offset at
#' infinity
#' \cr - \code{pars2[2] == 2.2} 1/n dependence in speciation rate
#' \cr - \code{pars2[2] == 2.3} exponential dependence in speciation rate with
#' parameter x (= exponent)
#' \cr - \code{pars2[2] == 3} linear dependence in
#' extinction rate
#' \cr - \code{pars2[2] == 4} exponential dependence in
#' extinction rate
#' \cr - \code{pars2[2] == 4.1} variant of exponential
#' dependence in extinction rate with offset at infinity
#' \cr - \code{pars2[2] ==
#' 4.2} 1/n dependence in extinction rate
#' \cr - \code{pars2[2] == 5} linear
#' dependence in speciation and extinction rate
#' \cr \cr \code{pars2[3]} sets
#' the conditioning:
#' \cr - \code{pars2[3] == 0} conditioning on stem or crown
#' age
#' \cr - \code{pars2[3] == 1} conditioning on stem or crown age and
#' non-extinction of the phylogeny
#' \cr - \code{pars2[3] == 2} conditioning on
#' stem or crown age and on the total number of extant taxa (including missing
#' species)
#' \cr - \code{pars2[3] == 3} conditioning on the total number of
#' extant taxa (including missing species)
#' \cr \cr \code{pars2[4]} sets whether
#' the likelihood is for the branching times (0) or the phylogeny (1)
#' \cr \cr \code{pars2[5]} sets whether the parameters and likelihood should be shown
#' on screen (1) or not (0)
#' \cr \cr \code{pars2[6]} sets whether the first data
#' point is stem age (1) or crown age (2)
#' @param brts A set of branching times of a phylogeny, all positive
#' @param missnumspec The number of species that are in the clade but missing
#' in the phylogeny
#' @param methode The method used to solve the master equation, default is
#' 'analytical' which uses matrix exponentiation; alternatively numerical ODE
#' solvers can be used. Before version 3.1 these were solvers from the deSolve
#' package such as 'lsoda' and 'ode45'. Currently solvers from odeint are used,
#' such as 'odeint::runge_kutta_cash_karp54', 'odeint::runge_kutta_fehlberg78',
#' 'odeint::runge_kutta_dopri5', or odeint::bulirsch_stoer'. The first two are
#' recommended in most cases.
#' @return The loglikelihood
#' @author Rampal S. Etienne & Bart Haegeman
#' @seealso \code{\link{dd_ML}}, \code{\link{dd_SR_loglik}},
#' \code{\link{dd_KI_loglik}}
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#' dd_loglik(pars1 = c(0.5,0.1,100), pars2 = c(100,1,1,1,0,2), brts = 1:10, missnumspec = 0) 
#' @export dd_loglik
dd_loglik = function(pars1,pars2,brts,missnumspec,methode = 'analytical')
{
   if(pars2[3] == 3)
   {
     rhs_func_name = 'dd_loglik_bw_rhs'
   } else
   {
     rhs_func_name = 'dd_loglik_rhs'
   }
   if(methode == 'analytical')
   {
     out = dd_loglik2(pars1,pars2,brts,missnumspec)
   } else
   {
     out = dd_loglik1(pars1,pars2,brts,missnumspec,methode = methode,rhs_func_name = rhs_func_name)
   }
   return(out)
}

dd_loglik1 = function(pars1,pars2,brts,missnumspec,methode = 'odeint::runge_kutta_cash_karp54',rhs_func_name = 'dd_loglik_rhs')
{
  if(length(pars2) == 4)
  {
    pars2[5] = 0
    pars2[6] = 2
  }
  ddep = pars2[2]
  cond = pars2[3]
  btorph = pars2[4]
  verbose = pars2[5]
  soc = pars2[6]
  if(cond == 3) { soc = 2 }
  la = pars1[1]
  mu = pars1[2]
  K = pars1[3]
  if(ddep == 5) {r = pars1[4]} else {r = 0}
  if(ddep == 1 | ddep == 5)
  { 
    lx = min(max(1 + missnumspec,1 + ceiling(la/(la - mu) * (r + 1) * K)),ceiling(pars2[1]))
  } else {
    if(ddep == 1.3)
    {
      lx = min(ceiling(K),ceiling(pars2[1]))
    } else {
      lx = round(pars2[1])
    }
  }
  if((ddep == 1) & ((mu == 0 & missnumspec == 0 & floor(K) != ceiling(K) & la > 0.05) | K == Inf))
  {
    loglik = bd_loglik(pars1[1:(2 + (K < Inf))],c(2*(mu == 0 & K < Inf),pars2[3:6]),brts,missnumspec)
  } else {
    abstol = 1e-10
    reltol = 1e-8 
    brts = -sort(abs(as.numeric(brts)),decreasing = TRUE)
    if(sum(brts == 0) == 0)
    {
      brts[length(brts) + 1] = 0
    }
    S = length(brts) + (soc - 2)
    if(min(pars1) < 0)
    {
      if(verbose) cat('The parameters are negative.\n')
      loglik = -Inf
    } else {
      if((mu == 0 & (ddep == 2 | ddep == 2.1 | ddep == 2.2)) | (la == 0 & (ddep == 4 | ddep == 4.1 | ddep == 4.2)) | (la <= mu))
      { 
        if(verbose) cat("These parameter values cannot satisfy lambda(N) = mu(N) for a positive and finite N.\n")
        loglik = -Inf
      } else {
        if(((ddep == 1 | ddep == 5) & ceiling(la/(la - mu) * (r + 1) * K) < (S + missnumspec)) | ((ddep == 1.3) & (S + missnumspec > ceiling(K))))
        {
          if(verbose) cat('The parameters are incompatible.\n')
          loglik = -Inf
        } else {
          loglik = (btorph == 0) * lgamma(S)
          if(cond != 3)
          {
            probs = rep(0,lx)
            probs[1] = 1 # change if other species at stem/crown age 
            for(k in 2:(S + 2 - soc))
            {
              k1 = k + (soc - 2)
              y = dd_integrate(probs,brts[(k-1):k],rhs_func_name,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
              probs = y[2,2:(lx+1)]
              if(is.na(sum(probs)) && pars1[2]/pars1[1] < 1E-4 && missnumspec == 0)
              { 
                loglik = dd_loglik_high_lambda(pars1 = pars1,pars2 = pars2,brts = brts)
                if(verbose) cat('High lambda approximation has been applied.\n')
                return(loglik)
              }
              if(k < (S + 2 - soc))
              {
                probs = flavec(ddep,la,mu,K,r,lx,k1) * probs # speciation event
              }
              cp <- check_probs(loglik,probs,verbose); loglik <- cp[[1]]; probs <- cp[[2]];
            }    
          } else {
            probs = rep(0,lx + 1)
            probs[1 + missnumspec] = 1
            for(k in (S + 2 - soc):2)
            {
              k1 = k + (soc - 2)
              y = dd_integrate(probs,-brts[k:(k-1)],rhs_func_name,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
              probs = y[2,2:(lx+2)]
              if(k > soc)
              {
                probs = c(flavec(ddep,la,mu,K,r,lx,k1-1),1) * probs # speciation event
              }
              cp <- check_probs(loglik,probs[1:lx],verbose); loglik <- cp[[1]]; probs[1:lx] <- cp[[2]];
            }
          }
          if(probs[1 + missnumspec] <= 0 | loglik == -Inf | is.na(loglik) | is.nan(loglik))
          {
            if(verbose) cat('Probabilities smaller than 0 or other numerical problems are encountered in final result.\n')
            loglik = -Inf
          } else {        
            loglik = loglik + (cond != 3 | soc == 1) * log(probs[1 + (cond != 3) * missnumspec]) - lgamma(S + missnumspec + 1) + lgamma(S + 1) + lgamma(missnumspec + 1)
            
            logliknorm = 0
            if(cond == 1 | cond == 2)
            {
              probsn = rep(0,lx)
              probsn[1] = 1 # change if other species at stem or crown age
              k = soc
              t1 = brts[1] 
              t2 = brts[S + 2 - soc]
              y = dd_integrate(probsn,c(t1,t2),rhs_func_name,c(pars1,k,ddep),rtol = reltol,atol = abstol,method = methode);
              probsn = y[2,2:(lx+1)]
              if(soc == 1) { aux = 1:lx }
              if(soc == 2) { aux = (2:(lx+1)) * (3:(lx+2))/6 }
              probsc = probsn/aux
              cp <- check_probs(logliknorm,probsc,verbose); logliknorm <- cp[[1]]; probsc <- cp[[2]];
              if(cond == 1) { logliknorm = logliknorm + log(sum(probsc)) }
              if(cond == 2) { logliknorm = logliknorm + log(probsc[S + missnumspec - soc + 1])}             
            }
            if(cond == 3)
            { 
              probsn = rep(0,lx + 1)
              probsn[S + missnumspec + 1] = 1
              
              TT = 1e14 # max(1,1/abs(la - mu)) * 1E+10 * max(abs(brts)) # make this more efficient later
              y = dd_integrate(probsn,c(0,TT),rhs_func_name,c(pars1,0,ddep),rtol = reltol,atol = abstol,method = methode)
              logliknorm = log(y[2,lx + 2])
              if(soc == 2)
              {
                probsn = rep(0,lx + 1)
                probsn[1:lx] = probs[1:lx]
                probsn = c(flavec(ddep,la,mu,K,r,lx,1),1) * probsn # speciation event
                y = dd_integrate(probsn,c(max(abs(brts)),TT),rhs_func_name,c(pars1,1,ddep),rtol = reltol,atol = abstol,method = methode)
                logliknorm = logliknorm - log(y[2,lx + 2])
              }
            }
            if(is.na(logliknorm) | is.nan(logliknorm) | logliknorm == Inf)
            {
              if(verbose) cat('The normalization did not yield a number.\n')
              loglik = -Inf
            } else
            {
              loglik = loglik - logliknorm
            }
          }
        }
      }
    }
    if(verbose)
    {
      s1 = sprintf('Parameters: %f %f %f',pars1[1],pars1[2],pars1[3])
      if(ddep == 5) {s1 = sprintf('%s %f',s1,pars1[4])}
      s2 = sprintf(', Loglikelihood: %f',loglik)
      cat(s1,s2,"\n",sep = "")
      utils::flush.console()
    }
  }
  loglik = as.numeric(loglik)
  if(is.nan(loglik) | is.na(loglik))
  {
    loglik = -Inf
  }
  return(loglik)
}

dd_loglik2 = function(pars1,pars2,brts,missnumspec)
{
if(length(pars2) == 4)
{
    pars2[5] = 0
    pars2[6] = 2
}
ddep = pars2[2]
cond = pars2[3]
btorph = pars2[4]
verbose <- pars2[5]
soc = pars2[6]
if(cond == 3)
{ 
    soc = 2
}
la = pars1[1]
mu = pars1[2]
K = pars1[3]
if(ddep == 5)
{
    r = pars1[4]
} else
{
    r = 0
}
if(ddep == 1 | ddep == 5)
{ 
    lx = min(max(1 + missnumspec,1 + ceiling(la/(la - mu) * (r + 1) * K)),ceiling(pars2[1]))
} else if(ddep == 1.3)
{
    lx = min(ceiling(K),ceiling(pars2[1]))
} else {
    lx = round(pars2[1])
}
if((ddep == 1) & ((mu == 0 & missnumspec == 0 & floor(K) != ceiling(K) & la > 0.05) | K == Inf))
{
    loglik = bd_loglik(pars1[1:(2 + (K < Inf))],c(2*(mu == 0 & K < Inf),pars2[3:6]),brts,missnumspec)
} else {
abstol = 1e-16
reltol = 1e-10 
brts = -sort(abs(as.numeric(brts)),decreasing = TRUE)
if(sum(brts == 0) == 0)
{
   brts[length(brts) + 1] = 0
}
S = length(brts) + (soc - 2)
if(min(pars1) < 0)
{
    loglik = -Inf
} else {
if((mu == 0 & (ddep == 2 | ddep == 2.1 | ddep == 2.2)) | (la == 0 & (ddep == 4 | ddep == 4.1 | ddep == 4.2)) | (la <= mu))
{ 
    if(verbose) cat("These parameter values cannot satisfy lambda(N) = mu(N) for a positive and finite N.\n")
    loglik = -Inf
} else {
    if(((ddep == 1 | ddep == 5) & ceiling(la/(la - mu) * (r + 1) * K) < (S + missnumspec)) | ((ddep == 1.3) & ((S + missnumspec) > ceiling(K))))
    {
       loglik = -Inf
    } else {
       loglik = (btorph == 0) * lgamma(S)
       if(cond != 3)
       {
          probs = rep(0,lx)
          probs[1] = 1 # change if other species at stem/crown age 
          for(k in 2:(S + 2 - soc))
          {
             k1 = k + (soc - 2)
             #y = deSolve::ode(probs,brts[(k-1):k],rhs_func,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
             #probs2 = y[2,2:(lx+1)]
             probs = dd_loglik_M(pars1,lx,k1,ddep,tt = abs(brts[k] - brts[k-1]),probs)
             if(is.na(sum(probs)) && pars1[2]/pars1[1] < 1E-4 && missnumspec == 0)
             {
               loglik = dd_loglik_high_lambda(pars1 = pars1,pars2 = pars2,brts = brts)
               if(verbose) cat('High lambda approximation has been applied.\n')
               return(loglik)
             }
             if(k < (S + 2 - soc))
             {
                 #probs = flavec(ddep,la,mu,K,r,lx,k1) * probs # speciation event
                 probs = lambdamu(0:(lx - 1) + k1,c(pars1[1:3],r),ddep)[[1]] * probs
             }
             cp <- check_probs(loglik,probs,verbose); loglik <- cp[[1]]; probs<- cp[[2]];
          }    
       } else {
          probs = rep(0,lx + 1)
          probs[1 + missnumspec] = 1
          for(k in (S + 2 - soc):2)
          {
             k1 = k + (soc - 2)
             #y = deSolve::ode(probs,-brts[k:(k-1)],dd_loglik_bw_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
             #probs2 = y[2,2:(lx+2)]
             probs = dd_loglik_M_bw(pars1,lx,k1,ddep,tt = abs(brts[k] - brts[k-1]),probs[1:lx])
             probs = c(probs,0)
             if(k > soc)
             {
                 #probs = c(flavec(ddep,la,mu,K,r,lx,k1-1),1) * probs # speciation event
                 probs = c(lambdamu(0:(lx - 1) + k1 - 1,pars1,ddep)[[1]],1) * probs
             }    
             cp <- check_probs(loglik,probs[1:lx],verbose); loglik <- cp[[1]]; probs[1:lx] <- cp[[2]];
          }
       }
       if(probs[1 + (cond != 3) * missnumspec] <= 0 | loglik == -Inf)
       {
          loglik = -Inf
       } else  {        
          loglik = loglik + (cond != 3 | soc == 1) * log(probs[1 + (cond != 3) * missnumspec]) - lgamma(S + missnumspec + 1) + lgamma(S + 1) + lgamma(missnumspec + 1)
  
          logliknorm = 0
          if(cond == 1 | cond == 2)
          {
             probsn = rep(0,lx)
             probsn[1] = 1 # change if other species at stem or crown age
             k = soc
             t1 = brts[1] 
             t2 = brts[S + 2 - soc]
             #y = deSolve::ode(probsn,c(t1,t2),rhs_func,c(pars1,k,ddep),rtol = reltol,atol = abstol,method = methode);
             #probsn = y[2,2:(lx+1)]
             probsn = dd_loglik_M(pars1,lx,k,ddep,tt = abs(t2 - t1),probsn)
             if(soc == 1) { aux = 1:lx }
             if(soc == 2) { aux = (2:(lx+1)) * (3:(lx+2))/6 }
             probsc = probsn/aux
             if(cond == 1) { logliknorm = log(sum(probsc)) }
             if(cond == 2) { logliknorm = log(probsc[S + missnumspec - soc + 1])}             
          }
          if(cond == 3)
          { 
             #probsn = rep(0,lx + 1)
             #probsn[S + missnumspec + 1] = 1 #/ (S + missnumspec)
             #TT = max(1,1/abs(la - mu)) * 100000000 * max(abs(brts)) # make this more efficient later
             #y = deSolve::ode(probsn,c(0,TT),dd_loglik_bw_rhs,c(pars1,0,ddep),rtol = reltol,atol = abstol,method = methode)
             #logliknorm = log(y[2,lx + 2])
             probsn = rep(0,lx + 1)
             probsn[2] = 1
             MM = dd_loglik_M_aux(pars1,lx + 1,k = 0,ddep)
             MM = MM[-1,-1]
             #probsn = SparseM::solve(-MM,probsn[2:(lx + 1)])
             MMinv = SparseM::solve(MM)
             probsn = -MMinv %*% probsn[2:(lx + 1)]
             logliknorm = log(probsn[S + missnumspec])
             if(soc == 2)
             {
                #probsn = rep(0,lx + 1)
                #probsn[1:lx] = probs[1:lx]
                #probsn = c(flavec(ddep,la,mu,K,r,lx,1),1) * probsn # speciation event
                #probsn = c(lambdamu(0:(lx - 1) + 1,pars1,ddep)[[1]],1) * probsn # speciation event
                #y = deSolve::ode(probsn,c(max(abs(brts)),TT),dd_loglik_bw_rhs,c(pars1,1,ddep),rtol = reltol,atol = abstol,method = methode)
                #logliknorm = logliknorm - log(y[2,lx + 2])
                probsn2 = rep(0,lx)
                probsn2 = lambdamu(0:(lx - 1) + 1,pars1,ddep)[[1]] * probs[1:lx]
                MM = dd_loglik_M_bw_aux(pars1,lx,k = 1,ddep)
                MMinv = SparseM::solve(MM)
                #probsn2 = SparseM::solve(-MM,probsn2[1:lx])
                probsn2 = -MMinv %*% probsn2[1:lx]
                logliknorm = logliknorm - log(probsn2[1])            
             }
          }
          loglik = loglik - logliknorm
       }
    }
}}
if(verbose)
{
    s1 = sprintf('Parameters: %f %f %f',pars1[1],pars1[2],pars1[3])
    if(ddep == 5) {s1 = sprintf('%s %f',s1,pars1[4])}
    s2 = sprintf(', Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    utils::flush.console()
}
}
loglik = as.numeric(loglik)
if(is.nan(loglik) | is.na(loglik) | loglik == Inf)
{
    loglik = -Inf
}
return(loglik)
}

dd_int <- function(initprobs,tvec,rhs_func,pars,rtol,atol,method)
{
  if(method == 'analytical')
  {
    lp <- length(pars)
    lx <- length(initprobs)
    ddep <- pars[lp]
    k1 <- pars[lp - 1]
    pars1 <- pars[-c(lp - 1,lp)]
    int <- dd_loglik_M(pars1,lx,k1,ddep,tvec[2] - tvec[1],initprobs)
  } else
  {
    int <- dd_integrate(initprobs,tvec,rhs_func,pars,rtol,atol,method)
  }
  return(int)
}

dd_integrate <- function(initprobs,tvec,rhs_func,pars,rtol,atol,method)
{
  if(method == 'analytical')
  {
    probs <- dd_loglik_M(pars = pars[1:3],
                         lx = length(initprobs),
                         k = pars[4],
                         ddep = pars[5],
                         tt = abs(tvec[2] - tvec[1]),
                         initprobs)
    y <- cbind(c(NA,NA),rbind(rep(NA,length(probs)),probs))
  } else
  {  
    #rhs_func_name <- 'no_name'
    #if(is.character(rhs_func))
    #{
      rhs_func_name <- rhs_func
      #if(rhs_func_name != 'dd_loglik_rhs' & rhs_func_name != 'dd_loglik_bw_rhs')
      #{
        rhs_func = match.fun(rhs_func)
      #}
    #}
    if(rhs_func_name == 'dd_loglik_rhs' || rhs_func_name == 'dd_loglik_bw_rhs')
    {
      parsvec = c(dd_loglik_rhs_precomp(pars,initprobs),pars[length(pars) - 1])
    } else 
    {
      parsvec = pars
    }
    if (startsWith(method, "odeint::")) {
      # couldn't find a better place to plug this in
      y <- dd_ode_odeint(initprobs, tvec, parsvec, atol, rtol, method = method, rhs_name = rhs_func_name)
    } else
    {
      y <- deSolve::ode(initprobs,tvec,rhs_func,parsvec,rtol = rtol,atol = atol,method = method)
    }
  }
  return(y)
}


dd_rhs_odeint_map = list(
  'dd_loglik_rhs' = dd_integrate_odeint,
  'dd_loglik_bw_rhs' = dd_integrate_bw_odeint
)


dd_ode_bw_right_open_odeint = function(initprobs, tvec, parsvec, method)
{
  y = initprobs
  dy = 0.0
  ddy = 0.0
  while (tvec[1] + 100 < tvec[2]) {
    y0 = y[length(y)]
    y <- dd_integrate_bw_odeint(y, c(tvec[1], tvec[1] + 100), parsvec, 1e-4, 1e-4, method)
    G = y[length(y)]
    dy = G - y0
    ddy = dy - ddy;
    tvec[1] = tvec[1] + 100
    if (abs(ddy) < 1e-10) {
      # dy became constant - extrapolate
      return(G + (tvec[2] - tvec[1]) * dy / 100.0)
    }
    if (dy > 1e20) {
      # blown up
      return(Inf)
    }
    ddy = dy
  }
  # never reached so far...
  y <- dd_integrate_bw_odeint(y, c(tvec[1], tvec[2]), parsvec, 1e-6, 1e-6, method)
  return(y[length(y)])
}


dd_ode_odeint = function(initprobs, tvec, parsvec, atol, rtol, method, rhs_name)
{
  fun = dd_rhs_odeint_map[[rhs_name]];
  # caller expect a matrix type
  lx = length(initprobs)
  y = matrix(nrow=2, ncol = lx + 1)
  if (tvec[2] > 10000 && rhs_name == 'dd_loglik_bw_rhs') {
    # tends to hang if not handled 
    y[2,2:(lx+1)] <- dd_ode_bw_right_open_odeint(initprobs, tvec, parsvec, method)
  }
  else {
    y[2,2:(lx+1)] <- fun(initprobs, tvec, parsvec, atol, rtol, method)
  }
  return(y)
}



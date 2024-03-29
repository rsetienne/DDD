# brts = branching times (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# - pars1[1] = la = (initial) speciation rate
# - pars1[2] = mu = extinction rate
# - pars1[3] = K = carrying capacity
# - pars1[4] = la2 = (initial) speciation rate
# - pars1[5] = mu2 = extinction rate
# - pars1[6] = K2 = carrying capacity
# - pars1[7] = tshift = time of shift
# - pars2[1] = lx = length of ODE variable x
# - pars2[2] = ddep = diversity-dependent model, mode of diversity-dependence
#  . ddep == 1 : linear dependence in speciation rate
#  . ddep == 1.3 : linear dependence in speciation rate with parameter K'
#  . ddep == 2 : exponential dependence in speciation rate
#  . ddep == 2.1: variant with offset at infinity
#  . ddep == 2.2: 1/n dependence in speciation rate
#  . ddep == 2.3: exponential dependence in speciation rate with parameter x
#  . ddep == 3 : linear dependence in extinction rate
#  . ddep == 4 : exponential dependence in extinction rate
#  . ddep == 4.1: variant with offset at infinity
#  . ddep == 4.2: 1/n dependence in speciation rate
# - pars2[3] = cond = conditioning
#  . cond == 0 : no conditioning
#  . cond == 1 : conditioning on non-extinction of the phylogeny
#  . cond == 2 : conditioning on non-extinction of the phylogeny and on the total number of extant taxa (including missing species)
#  . cond == 3 : conditioning on the total number of extant taxa (including missing species)
# - pars2[4] = btorph = likelihood of branching times (0) or phylogeny (1), differ by a factor (S - 1)! where S is the number of extant species
# - pars2[5] = parameters and likelihood should be printed (1) or not (0)
# - pars2[6] = likelihood is for a tree with crown age (2) or stem age (1)
# missnumspec = number of missing species    
# methode = the method used in the numerical solving of the set of the ode's



#' Loglikelihood for diversity-dependent diversification models with a shift in
#' the parameters at time t = tshift
#' 
#' This function computes loglikelihood of a diversity-dependent
#' diversification model for a given set of branching times and parameter
#' values where the parameters are allowed to shift at time t = tshift
#' 
#' 
#' @param pars1 Vector of parameters: \cr \cr \code{pars1[1]} corresponds to
#' lambda (speciation rate) before the shift \cr \code{pars1[2]} corresponds to
#' mu (extinction rate) before the shift \cr \code{pars1[3]} corresponds to K
#' (clade-level carrying capacity) before the shift \cr \code{pars1[4]}
#' corresponds to lambda (speciation rate) after the shift \cr \code{pars1[5]}
#' corresponds to mu (extinction rate) after the shift \cr \code{pars1[6]}
#' corresponds to K (clade-level carrying capacity) after the shift \cr
#' \code{pars1[7]} corresponds to tshift (the time of shift)
#' @param pars2 Vector of model settings: \cr \cr \code{pars2[1]} sets the
#' maximum number of species for which a probability must be computed.  This
#' must be larger than 1 + missnumspec + length(brts). \cr \cr \code{pars2[2]}
#' sets the model of diversity-dependence: \cr - \code{pars2[2] == 1} linear
#' dependence in speciation rate with parameter K (= diversity where speciation
#' = extinction)\cr - \code{pars2[2] == 1.3} linear dependence in speciation
#' rate with parameter K' (= diversity where speciation = 0)\cr -
#' \code{pars2[2] == 2} exponential dependence in speciation rate with
#' parameter K (= diversity where speciation = extinction)\cr - \code{pars2[2]
#' == 2.1} variant of exponential dependence in speciation rate with offset at
#' infinity\cr - \code{pars2[2] == 2.2} 1/n dependence in speciation rate\cr -
#' \code{pars2[2] == 2.3} exponential dependence in speciation rate with
#' parameter x (= exponent)\cr - \code{pars2[2] == 3} linear dependence in
#' extinction rate \cr - \code{pars2[2] == 4} exponential dependence in
#' extinction rate \cr - \code{pars2[2] == 4.1} variant of exponential
#' dependence in extinction rate with offset at infinity\cr - \code{pars2[2] ==
#' 4.2} 1/n dependence in extinction rate\cr\cr \code{pars2[3]} sets the
#' conditioning: \cr - \code{pars2[3] == 0} no conditioning \cr -
#' \code{pars2[3] == 1} conditioning on non-extinction of the phylogeny \cr -
#' \code{pars2[3] == 2} conditioning on non-extinction of the phylogeny and on
#' the total number of extant taxa (including missing species) \cr \cr
#' \code{pars2[4]} sets whether the likelihood is for the branching times (0)
#' or the phylogeny (1) \cr \cr \code{pars2[5]} sets whether the parameters and
#' likelihood should be shown on screen (1) or not (0) \cr \cr \code{pars2[6]}
#' sets whether the first data point is stem age (1) or crown age (2)
#' @param brts A set of branching times of a phylogeny, all positive
#' @param missnumspec The number of species that are in the clade but missing
#' in the phylogeny
#' @param methode The method used to solve the master equation, default is
#' 'analytical' which uses matrix exponentiation; alternatively numerical ODE
#' solvers can be used, such as 'odeint::runge_kutta_cash_karp54'. These were used in the
#' package before version 3.1.
#' @return The loglikelihood
#' @author Rampal S. Etienne & Bart Haegeman
#' @seealso \code{\link{dd_SR_ML}}, \code{\link{dd_loglik}},
#' \code{\link{dd_KI_loglik}}
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#' dd_SR_loglik(pars1 = c(0.2,0.1,50,0.2,0.1,70,5), pars2 = c(100,1,1,1,0,2),
#'    brts = 1:10, missnumspec = 0) 
#' @export dd_SR_loglik
dd_SR_loglik = function(pars1,pars2,brts,missnumspec,methode = 'analytical')
{
   if(methode == 'analytical')
   {
       out = dd_SR_loglik2(pars1,pars2,brts,missnumspec)
   } else {
       out = dd_SR_loglik1(pars1,pars2,brts,missnumspec,methode = methode)
   }
   return(out)
}

dd_SR_loglik1 = function(pars1,pars2,brts,missnumspec,methode = 'odeint::runge_kutta_cash_karp54')
{
if(length(pars2) == 4)
{
    pars2[5] = 0
    pars2[6] = 2
}
verbose <- pars2[5]
ddep = pars2[2]
abstol = 1e-16
reltol = 1e-10 
brts = -sort(abs(as.numeric(brts)),decreasing = TRUE)
if(sum(brts == 0) == 0) { brts[length(brts) + 1] = 0 }
soc = pars2[6]
S = length(brts) + (soc - 2)
if(min(pars1) < 0 || -pars1[7] <= min(brts))
{
    loglik = -Inf
} else {
if(((pars1[2] == 0 || pars1[4] == 0) && (ddep == 2 | ddep == 2.1 | ddep == 2.2)) || ((pars1[1] == 0 || pars1[3] == 0) && (ddep == 4 | ddep == 4.1 | ddep == 4.2)) || pars1[1] <= pars1[2] || pars1[4] <= pars1[5])
{ 
    cat("These parameter values cannot satisfy lambda(N) = mu(N) for some finite N.\n")
    loglik = -Inf
} else {
    la = pars1[1]
    mu = pars1[2]
    K = pars1[3]
    la2 = pars1[4]
    mu2 = pars1[5]
    K2 = pars1[6]
    tshift = -pars1[7]
    if(sum(abs(brts - tshift) < 1E-14) == 1) { tshift = tshift - 1E-8 }
    kshift = 1 + max(which(brts < tshift))
    if(ddep == 1) 
    { 
       lx = min(max(1 + missnumspec,1 + ceiling(max(la/(la - mu) * K,la2/(la2 - mu2) * K2))),ceiling(pars2[1]))
    } else {
       if(ddep == 1.3)
       {
          lx = min(max(ceiling(K),ceiling(K2)),ceiling(pars2[1]))
       } else {
          lx = round(pars2[1])
       }     
    }
    cond = pars2[3]
    btorph = pars2[4]

    if((ddep == 1 & (ceiling(la/(la - mu) * K) < kshift | ceiling(la2/(la2 - mu2) * K2) < (S + missnumspec))) | 
       (ddep == 1.3 & (ceiling(K) < kshift | ceiling(K2) < (S + missnumspec))))
    { 
       loglik = -Inf
    } else {
       loglik = (btorph == 0) * lgamma(S)
       if(cond != 3)
       {
           probs = rep(0,lx)
           probs[1] = 1 # change if other species at stem/crown age   
      
           if(kshift > 2)
           {
              for(k in 2:(kshift-1))
              {
                 k1 = k + (soc - 2)
                 y = dd_integrate(probs,brts[(k-1):k],'dd_loglik_rhs',c(pars1[1:3],k1,ddep),rtol = reltol,atol = abstol, method = methode)
                 probs = y[2,2:(lx+1)]
                 if(k < (S + 2 - soc))
                 {
                     probs = flavec(ddep,la,mu,K,0,lx,k1) * probs # speciation event
                 }
              }
              cp <- check_probs(loglik,probs,verbose); loglik <- cp[[1]]; probs <- cp[[2]];
           }   
           k = kshift
           k1 = k + (soc - 2)
           y = dd_integrate(probs,c(brts[k-1],tshift),'dd_loglik_rhs',c(pars1[1:3],k1,ddep),rtol = reltol,atol = abstol, method = methode)
           probs = y[2,2:(lx+1)]
           y = dd_integrate(probs,c(tshift,brts[k]),'dd_loglik_rhs',c(pars1[4:6],k1,ddep),rtol = reltol,atol = abstol, method = methode)
           probs = y[2,2:(lx+1)] 
           if(k < (S + 2 - soc))
           {
               probs = flavec(ddep,la2,mu2,K2,0,lx,k1) * probs # speciation event
           }
           cp <- check_probs(loglik,probs,verbose); loglik <- cp[[1]]; probs <- cp[[2]];
           if((kshift + 1) <= (S + 2 - soc))
           {
              for(k in (kshift + 1):(S + 2 - soc))
              {
                 k1 = k + (soc - 2)
                 y = dd_integrate(probs,brts[(k-1):k],'dd_loglik_rhs',c(pars1[4:6],k1,ddep),rtol = reltol,atol = abstol, method = methode)
                 probs = y[2,2:(lx+1)]
                 if(k < (S + 2 - soc))
                 {
                     probs = flavec(ddep,la2,mu2,K2,0,lx,k1) * probs # speciation event
                 }
                 cp <- check_probs(loglik,probs,verbose); loglik <- cp[[1]]; probs <- cp[[2]];
              }
           }    
       } else {
           probs = rep(0,lx + 1)
           probs[1 + missnumspec] = 1
           if((kshift + 1) <= S + 2 - soc)  
           {      
              for(k in (S + 2 - soc):(kshift + 1))
              {
                 k1 = k + (soc - 2)
                 y = deSolve::ode(probs,-brts[k:(k-1)],dd_loglik_bw_rhs,c(pars1[4:6],k1,ddep),rtol = reltol,atol = abstol, method = methode)
                 probs = y[2,2:(lx+2)]
                 if(k > soc)
                 {
                     probs = c(flavec(ddep,la2,mu2,K2,0,lx,k1-1),1) * probs # speciation event
                 }
                 cp <- check_probs(loglik,probs[1:lx],verbose); loglik <- cp[[1]]; probs[1:lx] <- cp[[2]];
              }
           }
           k = kshift
           k1 = k + (soc - 2)
           y = deSolve::ode(probs,-c(brts[k],tshift),dd_loglik_bw_rhs,c(pars1[4:6],k1,ddep),rtol = reltol,atol = abstol, method = methode)
           probs = y[2,2:(lx+2)]
           y = deSolve::ode(probs,-c(tshift,brts[k-1]),dd_loglik_bw_rhs,c(pars1[1:3],k1,ddep),rtol = reltol,atol = abstol, method = methode)
           probs = y[2,2:(lx+2)]
           if(k > soc)
           {
              probs = c(flavec(ddep,la,mu,K,0,lx,k1-1),1) * probs # speciation event
           }
           cp <- check_probs(loglik,probs[1:lx],verbose); loglik <- cp[[1]]; probs[1:lx] <- cp[[2]];
           for(k in (kshift-1):2)
           {
              k1 = k + (soc - 2)
              y = deSolve::ode(probs,-brts[k:(k-1)],dd_loglik_bw_rhs,c(pars1[1:3],k1,ddep),rtol = reltol,atol = abstol, method = methode)
              probs = y[2,2:(lx+2)]
              if(k > soc)
              {
                  probs = c(flavec(ddep,la,mu,K,0,lx,k1-1),1) * probs # speciation event
              }
              cp <- check_probs(loglik,probs[1:lx],verbose); loglik <- cp[[1]]; probs[1:lx] <- cp[[2]];
           }
       }
       if(probs[1 + (cond != 3) * missnumspec] <= 0 || loglik == -Inf)
       { 
          loglik = -Inf
       } else {
          loglik = loglik + (cond != 3 & soc == 2) * log(probs[1 + (cond != 3) * missnumspec]) - lgamma(S + missnumspec + 1) + lgamma(S + 1) + lgamma(missnumspec + 1)
                    
          logliknorm = 0
          if(cond == 1 | cond == 2)
          { 
             probs = rep(0,lx)
             probs[1] = 1 # change if other species at crown age
             k = soc
             y = dd_integrate(probs,c(brts[1],tshift),'dd_loglik_rhs',c(pars1[1:3],k,ddep),rtol = reltol,atol = abstol, method = methode);
             probs = y[2,2:(lx+1)]
             y = dd_integrate(probs,c(tshift,brts[length(brts)]),'dd_loglik_rhs',c(pars1[4:6],k,ddep),rtol = reltol,atol = abstol, method = methode);
             probs = y[2,2:(lx+1)]
             if(soc == 1) { aux = 1:lx }
             if(soc == 2) { aux = (2:(lx+1)) * (3:(lx+2))/6 }
             probsc = probs/aux
             if(cond == 1) { logliknorm = log(sum(probsc)) }
             if(cond == 2) { logliknorm = log(probsc[S + missnumspec - 1])}             
          }
          if(cond == 3)
          { 
             probsn = rep(0,lx + 1)
             probsn[S + missnumspec + 1] = 1
             TT = max(1,1/abs(la - mu),1/abs(la2 - mu2)) * 1E+8 * max(abs(brts)) # make this more efficient later
             y = deSolve::ode(probsn,c(0,-tshift),dd_loglik_bw_rhs,c(pars1[4:6],0,ddep),rtol = reltol,atol = abstol, method = methode)
             probsn = y[2,2:(lx+2)]
             y = deSolve::ode(probsn,c(-tshift,TT),dd_loglik_bw_rhs,c(pars1[1:3],0,ddep),rtol = reltol,atol = abstol, method = methode)
             logliknorm = log(y[2,lx + 2])
             if(soc == 2)
             {
                probsn = rep(0,lx + 1)
                probsn[1:lx] = probs[1:lx]
                probsn = c(flavec(ddep,la,mu,K,0,lx,1),1) * probsn # speciation event
                y = deSolve::ode(probsn,c(max(abs(brts)),TT),dd_loglik_bw_rhs,c(pars1[1:3],1,ddep),rtol = reltol,atol = abstol, method = methode)
                logliknorm = logliknorm - log(y[2,lx + 2])
             }
          }
          loglik = loglik - logliknorm   
       }
    }
}}
if(verbose)
{
    s1 = sprintf('Parameters: %f %f %f %f %f %f %f, ',pars1[1],pars1[2],pars1[3],pars1[4],pars1[5],pars1[6],pars1[7])
    s2 = sprintf('Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    utils::flush.console()
}
loglik = as.numeric(loglik)
if(is.nan(loglik) | is.na(loglik))
{
    loglik = -Inf
}
return(loglik)
}

dd_SR_loglik2 = function(pars1,pars2,brts,missnumspec,methode = 'odeint::runge_kutta_cash_karp54')
{
if(length(pars2) == 4)
{
    pars2[5] = 0
    pars2[6] = 2
}
verbose <- pars2[5]
ddep = pars2[2]
abstol = 1e-16
reltol = 1e-10 
brts = -sort(abs(as.numeric(brts)),decreasing = TRUE)
if(sum(brts == 0) == 0) { brts[length(brts) + 1] = 0 }
soc = pars2[6]
S = length(brts) + (soc - 2)
if(min(pars1) < 0 || -pars1[7] <= min(brts))
{
    loglik = -Inf
} else {
if(((pars1[2] == 0 || pars1[4] == 0) && (ddep == 2 | ddep == 2.1 | ddep == 2.2)) || ((pars1[1] == 0 || pars1[3] == 0) && (ddep == 4 | ddep == 4.1 | ddep == 4.2)) || pars1[1] <= pars1[2] || pars1[4] <= pars1[5])
{ 
    cat("These parameter values cannot satisfy lambda(N) = mu(N) for some finite N.\n")
    loglik = -Inf
} else {
    la = pars1[1]
    mu = pars1[2]
    K = pars1[3]
    la2 = pars1[4]
    mu2 = pars1[5]
    K2 = pars1[6]
    tshift = -pars1[7]
    if(sum(abs(brts - tshift) < 1E-14) == 1) { tshift = tshift - 1E-8 }
    kshift = 1 + max(which(brts < tshift))
    if(ddep == 1) 
    { 
       lx = min(max(1 + missnumspec,1 + ceiling(max(la/(la - mu) * K,la2/(la2 - mu2) * K2))),ceiling(pars2[1]))
    } else {
       if(ddep == 1.3)
       {
          lx = min(max(ceiling(K),ceiling(K2)),ceiling(pars2[1]))
       } else {
          lx = round(pars2[1])
       }     
    }
    cond = pars2[3]
    btorph = pars2[4]

    if((ddep == 1 & (ceiling(la/(la - mu) * K) < kshift | ceiling(la2/(la2 - mu2) * K2) < (S + missnumspec))) | 
       (ddep == 1.3 & (ceiling(K) < kshift | ceiling(K2) < (S + missnumspec))))
    { 
       loglik = -Inf
    } else {
       loglik = (btorph == 0) * lgamma(S)
       if(cond != 3)
       {
           probs = rep(0,lx)
           probs[1] = 1 # change if other species at stem/crown age   
      
           if(kshift > 2)
           {
              for(k in 2:(kshift-1))
              {
                 k1 = k + (soc - 2)
                 probs = dd_loglik_M(pars1[1:3],lx,k1,ddep,tt = abs(brts[k] - brts[k-1]),probs)
                 if(k < (S + 2 - soc))
                 {
                     #probs = flavec(ddep,la,mu,K,0,lx,k1) * probs # speciation event
                     probs = lambdamu(0:(lx - 1) + k1,c(pars1[1:3],0),ddep)[[1]] * probs
                 }
                 cp <- check_probs(loglik,probs,verbose); loglik <- cp[[1]]; probs <- cp[[2]];
              }
           }   
           k = kshift
           k1 = k + (soc - 2)
           probs = dd_loglik_M(pars1[1:3],lx,k1,ddep,tt = abs(tshift - brts[k-1]),probs)
           probs = dd_loglik_M(pars1[4:6],lx,k1,ddep,tt = abs(brts[k] - tshift),probs)
           if(k < (S + 2 - soc))
           {
               #probs = flavec(ddep,la2,mu2,K2,0,lx,k1) * probs # speciation event
               probs = lambdamu(0:(lx - 1) + k1,c(pars1[4:6],0),ddep)[[1]] * probs
           }
           cp <- check_probs(loglik,probs,verbose); loglik <- cp[[1]]; probs <- cp[[2]];
           if((kshift + 1) <= (S + 2 - soc))
           {
              for(k in (kshift + 1):(S + 2 - soc))
              {
                 k1 = k + (soc - 2)
                 probs = dd_loglik_M(pars1[4:6],lx,k1,ddep,tt = abs(brts[k] - brts[k - 1]),probs)
                 if(k < (S + 2 - soc))
                 {
                     #probs = flavec(ddep,la2,mu2,K2,0,lx,k1) * probs # speciation event
                     probs = lambdamu(0:(lx - 1) + k1,c(pars1[4:6],0),ddep)[[1]] * probs
                 }
                 cp <- check_probs(loglik,probs,verbose); loglik <- cp[[1]]; probs <- cp[[2]];
              }
           }    
       } else {
           probs = rep(0,lx + 1)
           probs[1 + missnumspec] = 1
           if((kshift + 1) <= S + 2 - soc)  
           {      
              for(k in (S + 2 - soc):(kshift + 1))
              {
                 k1 = k + (soc - 2)
                 #y = deSolve::ode(probs,-brts[k:(k-1)],dd_loglik_bw_rhs,c(pars1[4:6],k1,ddep),rtol = reltol,atol = abstol, method = methode)
                 #probs = y[2,2:(lx+2)]
                 probs = dd_loglik_M_bw(pars1[4:6],lx,k1,ddep,tt = abs(brts[k - 1] - brts[k]),probs[1:lx])
                 probs = c(probs,0) 
                 if(k > soc)
                 {
                     #probs = c(flavec(ddep,la2,mu2,K2,0,lx,k1-1),1) * probs # speciation event
                     probs = c(lambdamu(0:(lx - 1) + k1 - 1,pars1[4:6],ddep)[[1]],1) * probs
                 }
                 cp <- check_probs(loglik,probs[1:lx],verbose); loglik <- cp[[1]]; probs[1:lx] <- cp[[2]];
              }
           }
           k = kshift
           k1 = k + (soc - 2)
           #y = deSolve::ode(probs,-c(brts[k],tshift),dd_loglik_bw_rhs,c(pars1[4:6],k1,ddep),rtol = reltol,atol = abstol, method = methode)
           #probs = y[2,2:(lx+2)]
           probs = dd_loglik_M_bw(pars1[4:6],lx,k1,ddep,tt = abs(tshift - brts[k]),probs[1:lx])
           probs = c(probs,0)            
           #y = deSolve::ode(probs,-c(tshift,brts[k-1]),dd_loglik_bw_rhs,c(pars1[1:3],k1,ddep),rtol = reltol,atol = abstol, method = methode)
           #probs = y[2,2:(lx+2)]
           probs = dd_loglik_M_bw(pars1[1:3],lx,k1,ddep,tt = abs(brts[k - 1] - tshift),probs[1:lx])
           probs = c(probs,0) 
           if(k > soc)
           {
              #probs = c(flavec(ddep,la,mu,K,0,lx,k1-1),1) * probs # speciation event
              probs = c(lambdamu(0:(lx - 1) + k1 - 1,pars1[1:3],ddep)[[1]],1) * probs
           }
           cp <- check_probs(loglik,probs[1:lx],verbose); loglik <- cp[[1]]; probs[1:lx] <- cp[[2]];
           for(k in (kshift-1):2)
           {
              k1 = k + (soc - 2)
              #y = deSolve::ode(probs,-brts[k:(k-1)],dd_loglik_bw_rhs,c(pars1[1:3],k1,ddep),rtol = reltol,atol = abstol, method = methode)
              #probs = y[2,2:(lx+2)]
              probs = dd_loglik_M_bw(pars1[1:3],lx,k1,ddep,tt = abs(brts[k - 1] - brts[k]),probs[1:lx])
              probs = c(probs,0)
              if(k > soc)
              {
                  #probs = c(flavec(ddep,la,mu,K,0,lx,k1-1),1) * probs # speciation event
                  probs = c(lambdamu(0:(lx - 1) + k1 - 1,pars1[1:3],ddep)[[1]],1) * probs
              }
              cp <- check_probs(loglik,probs[1:lx],verbose); loglik <- cp[[1]]; probs[1:lx] <- cp[[2]];
           }
       }
       if(probs[1 + (cond != 3) * missnumspec] <= 0 || loglik == -Inf)
       { 
          loglik = -Inf
       } else {
          loglik = loglik + (cond != 3 & soc == 2) * log(probs[1 + (cond != 3) * missnumspec]) - lgamma(S + missnumspec + 1) + lgamma(S + 1) + lgamma(missnumspec + 1)
                    
          logliknorm = 0
          if(cond == 1 | cond == 2)
          { 
             probs = rep(0,lx)
             probs[1] = 1 # change if other species at crown age
             k = soc
             probs = dd_loglik_M(pars1[1:3],lx,k,ddep,tt = abs(tshift - brts[1]),probs)
             probs = dd_loglik_M(pars1[4:6],lx,k,ddep,tt = abs(brts[length(brts)] - tshift),probs)           
             if(soc == 1) { aux = 1:lx }
             if(soc == 2) { aux = (2:(lx+1)) * (3:(lx+2))/6 }
             probsc = probs/aux
             if(cond == 1) {
               sumprobsc <- sum(probsc)
               if(sumprobsc > 1) {
                 cat('Conditioning probability larger than 1 encountered; setting to 1.\n')
                 sumprobsc <- 1
               }
               logliknorm <- log(sum(probsc))
             }
             if(cond == 2) { logliknorm = log(probsc[S + missnumspec - 1])}             
          }
          if(cond == 3)
          { 
             #probsn = rep(0,lx + 1)
             #probsn[S + missnumspec + 1] = 1
             #TT = max(1,1/abs(la - mu),1/abs(la2 - mu2)) * 1E+10 * max(abs(brts)) # make this more efficient later
             #y = deSolve::ode(probsn,c(0,-tshift),dd_loglik_bw_rhs,c(pars1[4:6],0,ddep),rtol = reltol,atol = abstol, method = methode)
             #probsn = y[2,2:(lx+2)]
             #y = deSolve::ode(probsn,c(-tshift,TT),dd_loglik_bw_rhs,c(pars1[1:3],0,ddep),rtol = reltol,atol = abstol, method = methode)
             #logliknorm = log(y[2,lx + 2])
                       
             probsn = rep(0,lx + 1)
             probsn[S + missnumspec + 1] = 1
             MM2 = dd_loglik_M_bw_aux(pars1[4:6],lx + 1,k = 0,ddep)
             MM2inv = SparseM::solve(MM2[-1,-1])
             expMM2 = dd_loglik_M_bw(pars1[4:6],lx + 1,k = 0,ddep,tt = -tshift,probsn)
             MM1 = dd_loglik_M_bw_aux(pars1[1:3],lx + 1,k = 0,ddep)
             MM1inv = SparseM::solve(MM1[-1,-1])           
             expMM1 = dd_loglik_M_bw(pars1[1:3],lx + 1,k = 0,ddep,tt = -tshift,probsn)
             probsn = -MM2inv %*% probsn[2:(lx + 1)] + MM2inv %*% expMM2[2:(lx + 1)] - MM1inv %*% expMM2[2:(lx+ 1)]
             logliknorm = log(probsn[1])
             if(soc == 2)
             {
                #probsn = rep(0,lx + 1)
                #probsn[1:lx] = probs[1:lx]
                #probsn = c(flavec(ddep,la,mu,K,0,lx,1),1) * probsn # speciation event
                #y = deSolve::ode(probsn,c(max(abs(brts)),TT),dd_loglik_bw_rhs,c(pars1[1:3],1,ddep),rtol = reltol,atol = abstol, method = methode)
                #logliknorm = logliknorm - log(y[2,lx + 2])
                #print(log(y[2,lx + 2]))
                probsn2 = rep(0,lx)
                probsn2 = lambdamu(0:(lx - 1) + 1,pars1[1:3],ddep)[[1]] * probs[1:lx]
                MM = dd_loglik_M_bw_aux(pars1[1:3],lx,k = 1,ddep)
                MMinv = SparseM::solve(MM)
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
    s1 = sprintf('Parameters: %f %f %f %f %f %f %f, ',pars1[1],pars1[2],pars1[3],pars1[4],pars1[5],pars1[6],pars1[7])
    s2 = sprintf('Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    utils::flush.console()
}
loglik = as.numeric(loglik)
if(is.nan(loglik) | is.na(loglik))
{
    loglik = -Inf
}
return(loglik)
}

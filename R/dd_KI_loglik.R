# brtsM = branching times of main clade M (positive, from present to past)
# brtsS = branching times of subclade S (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# - pars1[1] = laM = (initial) speciation rate of main clade
# - pars1[2] = muM = extinction rate of main clade
# - pars1[3] = KM = carrying capacity of main clade
# - pars1[4] = laS = (initial) speciation rate of subclade
# - pars1[5] = muS = extinction rate of subclade
# - pars1[6] = KS = carrying capacity of subclade
# - pars1[7] = tinn = time of key innovation
# - pars2[1] = lx = length of ODE variable x
# - pars2[2] = ddep = diversity-dependent model, mode of diversity-dependence
#  . ddep == 1 : linear dependence in speciation rate with parameter K
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
# - pars2[4] = tsplit = time of split of innovative branch
# - pars2[5] = printing of parameters and likelihood (1) or not (0)
# - pars2[6] = likelihood is for a tree with crown age (2) or stem age (1)
# missnumspec = number of missing species in main clade M and subclade S
# methode = the method used in the numerical solving of the set of the ode's



#' Loglikelihood for diversity-dependent diversification models with decoupling
#' of a subclade from a main clade at time t = t_d
#' 
#' This function computes loglikelihood of a diversity-dependent
#' diversification model for a given set of branching times and parameter
#' values where the diversity-dependent dynamics of a subclade decouple from
#' the dynamics of the main clade at time t_d, potentially accompanied by a
#' shift in parameters.
#' 
#' 
#' @param pars1 Vector of parameters: \cr \cr \code{pars1[1]} corresponds to
#' lambda_M (speciation rate) of the main clade \cr \code{pars1[2]} corresponds
#' to mu_M (extinction rate) of the main clade \cr \code{pars1[3]} corresponds
#' to K_M (clade-level carrying capacity) of the main clade \cr \code{pars1[4]}
#' corresponds to lambda_M (speciation rate) of the subclade \cr
#' \code{pars1[5]} corresponds to mu_S (extinction rate) of the subclade \cr
#' \code{pars1[6]} corresponds to K_S (clade-level carrying capacity) of the
#' subclade \cr \code{pars1[7]} corresponds to t_d (the time of decoupling)
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
#' \code{pars2[3] == 1} conditioning on non-extinction of the phylogeny \cr \cr
#' \code{pars2[4]} sets the time of splitting of the branch that will decouple
#' \cr \cr \code{pars2[5]} sets whether the parameters and likelihood should be
#' shown on screen (1) or not (0) \cr \cr \code{pars2[6]} sets whether the
#' first data point is stem age (1) or crown age (2) \cr \cr \code{pars2[7]}
#' sets whether the old (incorrect) likelihood should be used (0), or whether
#' new corrected version should be used (1)
#' @param brtsM A set of branching times of the main clade in the phylogeny,
#' all positive
#' @param brtsS A set of branching times of the subclade in the phylogeny, all
#' positive
#' @param missnumspec The number of species that are in the clade but missing
#' in the phylogeny. One can specify the sum of the missing species in main
#' clade and subclade or a vector c(missnumspec_M,missnumspec_S) with missing
#' species in main clade and subclade respectively.
#' @param methode The method used to solve the master equation, default is
#' 'analytical' which uses matrix exponentiation; alternatively numerical ODE
#' solvers can be used, such as 'lsoda' or 'ode45'. These were used in the
#' package before version 3.1.
#' @return The loglikelihood
#' @author Rampal S. Etienne & Bart Haegeman
#' @seealso \code{\link{dd_KI_ML}}, \code{\link{dd_loglik}}
#' \code{\link{dd_SR_loglik}}
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#' 
#' pars1 = c(0.25,0.12,25.51,1.0,0.16,8.61,9.8)
#' pars2 = c(200,1,0,18.8,1,2)
#' missnumspec = 0
#' brtsM = c(25.2,24.6,24.0,22.5,21.7,20.4,19.9,19.7,18.8,17.1,15.8,11.8,9.7,8.9,5.7,5.2)
#' brtsS = c(9.6,8.6,7.4,4.9,2.5)
#' dd_KI_loglik(pars1,pars2,brtsM,brtsS,missnumspec,methode = 'ode45')
#' 
#' @export dd_KI_loglik
dd_KI_loglik = function(pars1,pars2,brtsM,brtsS,missnumspec,methode = 'analytical')
{ 
   if(length(pars2) == 4)
   {
      pars2[5] = 0
      pars2[6] = 2
      pars2[7] = 1
   }
   if(is.na(pars2[7]))
   {
     pars2[7] = 0
   }
   if(methode == 'analytical')
   {
       out = dd_KI_loglik2(pars1,pars2,brtsM,brtsS,missnumspec)
   } else {
       out = dd_KI_loglik1(pars1,pars2,brtsM,brtsS,missnumspec,methode = methode)
   }
   return(out)
}

dd_KI_loglik1 = function(pars1,pars2,brtsM,brtsS,missnumspec,methode = 'lsoda')
{
abstol = 1e-16
reltol = 1e-14
m = missnumspec
# order branching times
brts = -sort(abs(c(brtsM,brtsS)),decreasing = TRUE)
if(sum(brts == 0) == 0)
{ 
   brts[length(brts) + 1] = 0
}
verbose <- pars2[5]
soc = pars2[6]
S = length(brts) + (soc - 2)
brtsM = -sort(abs(brtsM),decreasing = TRUE)
if(sum(brtsM == 0) == 0)
{ 
   brtsM[length(brtsM) + 1] = 0
}
brtsS = -sort(abs(brtsS),decreasing = TRUE)
if(sum(brtsS == 0) == 0)
{
   brtsS[length(brtsS) + 1] = 0
}

if(min(pars1) < 0 | -pars1[7] <= min(brtsM) | -pars1[7] >= min(brtsS))
{
    loglik = -Inf
} else {
if(((pars1[2] == 0 || pars1[4] == 0) && pars2[2] == 2) | ((pars1[1] == 0 | pars1[3] == 0) & pars2[2] == 4) | pars1[1] <= pars1[2] | pars1[4] <= pars1[5])
{ 
    cat("These parameter values cannot satisfy lambda(N) = mu(N) for some finite N.\n")
    loglik = -Inf
} else {
    laM = pars1[1]
    muM = pars1[2]
    KM = pars1[3]
    laS = pars1[4]
    muS = pars1[5]
    KS = pars1[6]
    tinn = -pars1[7]
    lmax = pars2[1]
    ddep = pars2[2]
    if(ddep == 1)
    {
        lxM = min(max(1 + m[1],1 + ceiling(laM/(laM - muM) * KM)),ceiling(lmax))
        lxS = min(max(1 + m[1],1 + ceiling(laS/(laS - muS) * KS)),ceiling(lmax))
    } else {
       if(ddep == 1.3)
       {
          lxM = min(max(1 + m[1],1 + ceiling(KM)),ceiling(lmax))
          lxS = min(max(1 + m[1],1 + ceiling(KS)),ceiling(lmax))         
       } else {
          lxM = round(lmax)
          lxS = round(lmax)
       }
    }

    n0 = (ddep == 2 | ddep == 4)
    cond = pars2[3]
    tsplit = -pars2[4]    
    S1 = length(brtsM) - 1 + (soc - 2)
    if(sum(brtsS == tinn) == 0) { brtsS = c(tinn,brtsS) }
    S2 = length(brtsS) - 1
    S1a = S1
    S2a = S2
    summ = sum(m)
    if(length(m) == 2)
    {
       S1a = S1 + m[1]
       S2a = S2 + m[2]
       summ = 0
    }  
    if((ddep == 1 & ( (ceiling(laM/(laM - muM) * KM) < S1a) | (ceiling(laS/(laS - muS) * KS) < S2a) )) |
       (ddep == 1.3 & ( (ceiling(KM) < S1a) | (ceiling(KS) < S2a) | (ceiling(KM) + ceiling(KS) < S1a + S2a + summ) )))
    { 
       loglik = -Inf
    } else {

    # avoid coincidence of branching time and key innovation time
    if(sum(abs(brtsM - tinn) < 1E-14) == 1) { tinn = tinn - 1E-8 }

    # compute likelihood of clade M
    loglikM = 0
    lx = lxM
    probs = rep(0,lx)
    probs[1] = 1 # change if other species at crown age

    ka = sum(brtsM < tinn);
    for(k in 2:(ka+1))
    {
       k1 = k + (soc - 2)
       t1 = brtsM[k - 1]; t2 = min(c(tinn,brtsM[k]))
       y = dd_integrate(probs,c(t1,t2),'dd_loglik_rhs',c(pars1[1:3],k1,ddep),rtol = reltol,atol = abstol,method = methode)
       probs = y[2,2:(lx + 1)]
       if(t2 < tinn)
       {
           probs = flavec(ddep,laM,muM,KM,0,lxM,k1,n0) * probs # speciation event
       }
       cp <- check_probs(loglikM,probs,verbose); loglikM <- cp[[1]]; probs <- cp[[2]];
    }
    for(k in (ka + 1):max(ka + 1,S1 + 1))
    {
       k1 = k + (soc - 2)
       t1 = max(tinn,brtsM[k - 1]); t2 = brtsM[k];
       if(pars2[7] == 0) { probs = probs * k1/(k1 + (0:(length(probs) - 1))) }
       y = dd_integrate(probs,c(t1,t2),'dd_loglik_rhs',c(pars1[1:3],k1-1,ddep),rtol = reltol,atol = abstol,method = methode)
       probs = y[2,2:(lx + 1)]
       if(k < (S1+1))
       {
           probs = flavec(ddep,laM,muM,KM,0,lxM,k1-1,n0) * probs # speciation event
       }
       cp <- check_probs(loglikM,probs,verbose); loglikM <- cp[[1]]; probs <- cp[[2]];
    }
    if(length(m) == 1)
    { 
       loglikM = loglikM + log(probs[1 + (0:m)])   
    } else {
       loglikM = loglikM + log(probs[1 + m[1]])   
    }
    # compute likelihood of clade S
    loglikS = 0
    lx = lxS
    probs = rep(0,lx)
    probs[1] = 1
    for(k in 1:S2)
    {
       t1 = brtsS[k]; t2 = brtsS[k+1]
       y = dd_integrate(probs,c(t1,t2),'dd_loglik_rhs',c(pars1[4:6],k,ddep),rtol = reltol,atol = abstol,method = methode)
       probs = y[2,2:(lx+1)]
       if(k < S2)
       {
           probs = flavec(ddep,laS,muS,KS,0,lxS,k,n0) * probs # speciation event
       }
       cp <- check_probs(loglikS,probs,verbose); loglikS <- cp[[1]]; probs <- cp[[2]];
    }
    if(length(m) == 1)
    {
       loglikS = loglikS + log(probs[1 + (0:m)])
    } else {
       loglikS = loglikS + log(probs[1 + m[2]])
    }

    #if(cond == 3)
    #{
    #   loglikS = 0
    #   lx = lxS
    #   probs = rep(0,lx + 1)
    #   if(length(m) == 1)
    #   {
    #      probs[1:(m+1)] = 1
    #   } else {
    #      probs[1 + m[2]] = 1
    #   }
    #   for(k in (S2 + 1):2)
    #   {
    #      k1 = k - 1
    #      y = deSolve::ode(probs,-brts[k:(k-1)],dd_loglik_bw_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
    #      probs = y[2,2:(lx+2)]
    #      if(k1 > 1)
    #      {
    #         probs = c(flavec(ddep,la,mu,K,r,lx,k1-1,n0),1) * probs # speciation event
    #      }
    #      if(sum(probs[1:lx]) <= 0)
    #      {
    #         loglik = -Inf
    #         break
    #      } else {
    #         loglikS = loglikS + log(sum(probs[1:lx]))
    #      }
    #      probs[1:lx] = probs[1:lx]/sum(probs[1:lx])
    #   }    
    #   loglikS = loglikS + log(probs[1])
    #}

    # total likelihood = likelihood clade M x likelihood clade S
    if(length(m) == 1)
    {
       #loglik = log(sum(exp(loglikM + loglikS[length(loglikS):1])))
       loglikMmax = max(loglikM)
       loglikSmax = max(loglikS)
       loglikMdelta = loglikM - loglikMmax
       loglikSdelta = loglikS - loglikSmax
       loglik = loglikMmax + loglikSmax + log(sum(exp(loglikMdelta + loglikSdelta[length(loglikSdelta):1])))
    } else {
       loglik = loglikM + loglikS
    }
    if(is.nan(loglik) | is.na(loglik))
    {
       loglik = -Inf
    }   

    if(cond == 0 | loglik == -Inf)
    {
       logliknorm = 0
    } else {   
       # COMPUTE NORMALIZATION
       tcrown = brts[1]
       tpres = 0
       # compute survival probability of clade S
       lx = lxS
       nx = -1:lx
       if(ddep == 1) 
       { 
           lavec = pmax(rep(0,lx + 2),laS - (laS - muS)/KS * nx)
           muvec = muS * rep(1,lx + 2)
       } 
       if(ddep == 1.3) 
       { 
           lavec = pmax(rep(0,lx + 2),laS * (1 - nx/KS))
           muvec = muS * rep(1,lx + 2)
       } 
       if(ddep == 2 | ddep == 2.1 | ddep == 2.2)
       {
           x = -(log(laS/muS)/log(KS+n0))^(ddep != 2.2) 
           lavec = pmax(rep(0,lx + 2),laS * (nx + n0)^x)
           muvec = muS * rep(1,lx + 2)
       }
       if(ddep == 2.3)
       {
           x = KS 
           lavec = pmax(rep(0,lx + 2),laS * (nx + n0)^x)
           muvec = muS * rep(1,lx + 2)
       }       
       if(ddep == 3)
       {
           lavec = laS * rep(1,lx + 2)
           muvec = muS + (laS - muS)/KS * nx
       }    
       if(ddep == 4 | ddep == 4.1 | ddep == 4.2)
       {
           lavec = laS * rep(1,lx + 2)
           x = (log(laS/muS)/log(KS+n0))^(ddep != 4.2)
           muvec = (nx + n0)^x
       }    
       m1 = lavec[1:lx] * nx[1:lx]
       m2 = muvec[3:(lx + 2)] * nx[3:(lx + 2)]
       m3 = (lavec[2:(lx + 1)] + muvec[2:(lx + 1)]) * nx[2:(lx + 1)]
       probs = rep(0,lx) # probs[1] = extinction probability
       probs[2] = 1 # clade S starts with one species
       y = deSolve::ode(probs,c(tinn,tpres),dd_logliknorm_rhs1,c(m1,m2,m3),rtol = reltol,atol = abstol,method = methode)
       probs = y[2,2:(lx+1)]   
       PS = 1 - probs[1]
   
       # compute survival probability of clade M
       lx = lxM
       nx1 = rep(-1:lx,lx + 2)
       dim(nx1) = c(lx + 2,lx + 2) # row index = number of species in first group 
       nx2 = t(nx1) # column index = number of species in second group
       nxt = nx1 + nx2
       if(ddep == 1) 
       { 
           lavec = pmax(matrix(0,lx + 2,lx + 2),laM - (laM-muM)/KM * nxt)
           muvec = muM * matrix(1,lx + 2,lx + 2)
       }         
       if(ddep == 1.3) 
       { 
           lavec = pmax(matrix(0,lx + 2,lx + 2),laM * (1 - nxt/KM))
           muvec = muM * matrix(1,lx + 2,lx + 2)
       } 
       if(ddep == 2 | ddep == 2.1 | ddep == 2.2)
       { 
           x = -(log(laM/muM)/log(KM+n0))^(ddep != 2.2)
           lavec = pmax(matrix(0,lx + 2,lx + 2),laM * (nxt + n0)^x)
           muvec = muM * matrix(1,lx + 2,lx + 2)
       }
       if(ddep == 2.3)
       { 
           x = KM
           lavec = pmax(matrix(0,lx + 2,lx + 2),laM * (nxt + n0)^x)
           muvec = muM * matrix(1,lx + 2,lx + 2)
       }    
       if(ddep == 3)
       {
           lavec = laM * matrix(1,lx + 2,lx + 2)
           muvec = muM + (laM - muM)/KM * nxt
       }    
       if(ddep == 4 | ddep == 4.1 | ddep == 4.2)
       {
           lavec = laM * matrix(1,lx + 2,lx + 2)
           x = (log(laM/muM)/log(KM+n0))^(ddep != 4.2)
           muvec = (nxt + n0)^x
       }    
       m1 = lavec[1:lx,2:(lx+1)] * nx1[1:lx,2:(lx+1)]
       m2 = muvec[3:(lx+2),2:(lx+1)] * nx1[3:(lx+2),2:(lx+1)]
       ma = lavec[2:(lx+1),2:(lx+1)] + muvec[2:(lx+1),2:(lx+1)]
       m3 = ma * nx1[2:(lx+1),2:(lx+1)]
       m4 = lavec[2:(lx+1),1:lx] * nx2[2:(lx+1),1:lx]
       m5 = muvec[2:(lx+1),3:(lx+2)] * nx2[2:(lx+1),3:(lx+2)]
       m6 = ma * nx2[2:(lx+1),2:(lx+1)]
       probs = matrix(0,lx,lx)

       # probs[1,1] = probability of extinction of both lineages
       # sum(probs[1:lx,1]) = probability of extinction of second lineage
       probs[2,2] = 1 # clade M starts with two species
       # STEP 1: integrate from tcrown to tinn
       dim(probs) = c(lx*lx,1)
       y = deSolve::ode(probs,c(tcrown,tinn),dd_logliknorm_rhs2,list(m1,m2,m3,m4,m5,m6),rtol = reltol,atol = abstol, method = "ode45")
       probs = y[2,2:(lx * lx + 1)]
       dim(probs) = c(lx,lx)
       probs[1,1:lx] = 0
       probs[1:lx,1] = 0
       # STEP 2: transformation at tinn
       nx1a = nx1[2:(lx+1),2:(lx+1)]
       nx2a = nx2[2:(lx+1),2:(lx+1)]
       probs = probs * nx1a/(nx1a+nx2a)
       probs = rbind(probs[2:lx,1:lx], rep(0,lx))
       dim(probs) = c(lx * lx,1)
       # STEP 3: integrate from tinn to tpres
       y = deSolve::ode(probs,c(tinn,tpres),dd_logliknorm_rhs2,list(m1,m2,m3,m4,m5,m6),rtol = reltol,atol = abstol, method = "ode45")
       probs = y[2,2:(lx * lx + 1)]
       dim(probs) = c(lx,lx)
       PM12 = sum(probs[2:lx,2:lx])
       PM2 = sum(probs[1,2:lx])
       logliknorm = log(2) + (cond == 1) * log(PM12 + PS * PM2) + (cond == 4) * (log(PM12 + PM2) + log(PS))
    }
    if(length(m) > 1)
    {
       Sv = c(S1,S2)
    } else {
       Sv = S
    }
    loglik = loglik - logliknorm - sum(lgamma(Sv + m + 1) - lgamma(Sv + 1) - lgamma(m + 1))
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


dd_KI_loglik2 = function(pars1,pars2,brtsM,brtsS,missnumspec,methode = 'lsoda')
{
abstol = 1e-16
reltol = 1e-14
m = missnumspec
# order branching times
brts = -sort(abs(c(brtsM,brtsS)),decreasing = TRUE)
if(sum(brts == 0) == 0)
{ 
   brts[length(brts) + 1] = 0
}
verbose <- pars2[5]
soc = pars2[6]
S = length(brts) + (soc - 2)
brtsM = -sort(abs(brtsM),decreasing = TRUE)
if(sum(brtsM == 0) == 0)
{ 
   brtsM[length(brtsM) + 1] = 0
}
brtsS = -sort(abs(brtsS),decreasing = TRUE)
if(sum(brtsS == 0) == 0)
{
   brtsS[length(brtsS) + 1] = 0
}

if(min(pars1) < 0 | -pars1[7] <= min(brtsM) | -pars1[7] >= min(brtsS))
{
    loglik = -Inf
} else {
if(((pars1[2] == 0 || pars1[4] == 0) && pars2[2] == 2) | ((pars1[1] == 0 | pars1[3] == 0) & pars2[2] == 4) | pars1[1] <= pars1[2] | pars1[4] <= pars1[5])
{ 
    cat("These parameter values cannot satisfy lambda(N) = mu(N) for some finite N.\n")
    loglik = -Inf
} else {
    laM = pars1[1]
    muM = pars1[2]
    KM = pars1[3]
    laS = pars1[4]
    muS = pars1[5]
    KS = pars1[6]
    tinn = -pars1[7]
    lmax = pars2[1]
    ddep = pars2[2]
    if(ddep == 1)
    {
        lxM = min(max(1 + m[1],1 + ceiling(laM/(laM - muM) * KM)),ceiling(lmax))
        lxS = min(max(1 + m[1],1 + ceiling(laS/(laS - muS) * KS)),ceiling(lmax))
    } else if(ddep == 1.3)
    {
        lxM = min(max(1 + m[1],1 + ceiling(KM)),ceiling(lmax))
        lxS = min(max(1 + m[1],1 + ceiling(KS)),ceiling(lmax))         
    } else {
        lxM = round(lmax)
        lxS = round(lmax)
    }
    n0 = (ddep == 2 | ddep == 4)
    cond = pars2[3]
    tsplit = -pars2[4]    
    S1 = length(brtsM) - 1 + (soc - 2)
    if(sum(brtsS == tinn) == 0)
    {
       brtsS = c(tinn,brtsS)
    }
    S2 = length(brtsS) - 1
    S1a = S1
    S2a = S2
    summ = sum(m)
    if(length(m) == 2)
    {
       S1a = S1 + m[1]
       S2a = S2 + m[2]
       summ = 0
    }  
    if((ddep == 1 & ( (ceiling(laM/(laM - muM) * KM) < S1a) | (ceiling(laS/(laS - muS) * KS) < S2a) )) |
       (ddep == 1.3 & ( (ceiling(KM) < S1a) | (ceiling(KS) < S2a) | (ceiling(KM) + ceiling(KS) < S1a + S2a + summ) )))
    { 
       loglik = -Inf
    } else {

    # avoid coincidence of branching time and key innovation time
    if(sum(abs(brtsM - tinn) < 1E-14) == 1)
    {
        tinn = tinn - 1E-8
    }

    # compute likelihood of clade M
    loglikM = 0
    lx = lxM
    probs = rep(0,lx)
    probs[1] = 1 # change if other species at crown age

    ka = sum(brtsM < tinn);
    for(k in 2:(ka+1))
    {
       k1 = k + (soc - 2)
       t1 = brtsM[k - 1]; t2 = min(c(tinn,brtsM[k]))
       #y = deSolve::ode(probs,c(t1,t2),dd_loglik_rhs,c(pars1[1:3],k1,ddep),rtol = reltol,atol = abstol,method = methode)
       #probs2 = y[2,2:(lx + 1)]
       #MM = dd_loglik_M_aux(pars1[1:3],lx,k1,ddep)
       #probs3 = Matrix::expm(MM * abs(t2 - t1)) %*% probs
       probs = dd_loglik_M(pars1[1:3],lx,k1,ddep,tt = abs(t2 - t1),probs)
       if(t2 < tinn)
       {
           probs = lambdamu(0:(lx - 1) + k1,c(pars1[1:3],0),ddep)[[1]] * probs
       }
       cp <- check_probs(loglikM,probs,verbose); loglikM <- cp[[1]]; probs <- cp[[2]];
    }
    for(k in (ka + 1):max(ka + 1,S1 + 1))
    {
       k1 = k + (soc - 2)
       t1 = max(tinn,brtsM[k - 1]); t2 = brtsM[k];
       if(pars2[7] == 1) { probs = probs * k1/(k1 + (0:(length(probs) - 1))) }
       #y = deSolve::ode(probs,c(t1,t2),dd_loglik_rhs,c(pars1[1:3],k1-1,ddep),rtol = reltol,atol = abstol,method = methode)
       #probs = y[2,2:(lx + 1)]
       probs = dd_loglik_M(pars1[1:3],lx,k1-1,ddep,tt = abs(t2 - t1),probs)
       if(k < (S1+1))
       {
           #probs = flavec(ddep,laM,muM,KM,0,lxM,k1-1,n0) * probs # speciation event
           probs = lambdamu(0:(lx - 1) + k1 - 1,c(pars1[1:3],0),ddep)[[1]] * probs
       }
       cp <- check_probs(loglikM,probs,verbose); loglikM <- cp[[1]]; probs <- cp[[2]];
    }
    if(length(m) == 1)
    { 
       loglikM = loglikM + log(probs[1 + (0:m)])   
    } else {
       loglikM = loglikM + log(probs[1 + m[1]])   
    }
    # compute likelihood of clade S
    loglikS = 0
    lx = lxS
    probs = rep(0,lx)
    probs[1] = 1
    for(k in 1:S2)
    {
       t1 = brtsS[k]; t2 = brtsS[k+1]
       #y = deSolve::ode(probs,c(t1,t2),dd_loglik_rhs,c(pars1[4:6],k,ddep),rtol = reltol,atol = abstol,method = methode)
       #probs = y[2,2:(lx+1)]
       probs = dd_loglik_M(pars1[4:6],lx,k,ddep,tt = abs(t2 - t1),probs)
       if(k < S2)
       {
           #probs = flavec(ddep,laS,muS,KS,0,lxS,k,n0) * probs # speciation event
           probs = lambdamu(0:(lx - 1) + k,c(pars1[4:6],0),ddep)[[1]] * probs
       }
       cp <- check_probs(loglikS,probs,verbose); loglikS <- cp[[1]]; probs <- cp[[2]];
    }
    if(length(m) == 1)
    {
       loglikS = loglikS + log(probs[1 + (0:m)])
    } else {
       loglikS = loglikS + log(probs[1 + m[2]])
    }

    # total likelihood = likelihood clade M x likelihood clade S
    if(length(m) == 1)
    {
       #loglik = log(sum(exp(loglikM + loglikS[length(loglikS):1])))
       loglikMmax = max(loglikM)
       loglikSmax = max(loglikS)
       loglikMdelta = loglikM - loglikMmax
       loglikSdelta = loglikS - loglikSmax
       loglik = loglikMmax + loglikSmax + log(sum(exp(loglikMdelta + loglikSdelta[length(loglikSdelta):1])))
    } else {
       loglik = loglikM + loglikS
    }

    if(is.nan(loglik) | is.na(loglik))
    {
       loglik = -Inf
    }   

    if(cond == 0 | loglik == -Inf)
    {
       logliknorm = 0
    } else {   
       # COMPUTE NORMALIZATION
       tcrown = brts[1]
       tpres = 0
       # compute survival probability of clade S
       lx = lxS
       probs = rep(0,lx) # probs[1] = extinction probability
       probs[2] = 1 # clade S starts with one species
       #nx = -1:lx
       #lambdamu = lambdamu(nx,c(pars1[4:6],0),ddep,tt = 1)
       #lavec = lambdamu[[1]]
       #muvec = lambdamu[[2]]
       #m1 = lavec[1:lx] * nx[1:lx]
       #m2 = muvec[3:(lx + 2)] * nx[3:(lx + 2)]
       #m3 = (lavec[2:(lx + 1)] + muvec[2:(lx + 1)]) * nx[2:(lx + 1)]
       #y = deSolve::ode(probs,c(tinn,tpres),dd_logliknorm_rhs1,c(m1,m2,m3),rtol = reltol,atol = abstol,method = methode)
       #probs = y[2,2:(lx+1)]   
       #PS = 1 - probs[1]
       probs = dd_loglik_M(pars1[4:6],lx,0,ddep,tt = abs(tpres - tinn),probs)
       PS = 1 - probs[1]
   
       # compute survival probability of clade M
       lx = lxM
       #nx1 = rep(-1:lx,lx + 2)
       #dim(nx1) = c(lx + 2,lx + 2) # row index = number of species in first group 
       #nx2 = t(nx1) # column index = number of species in second group
       #lambdamu = lambdamu2(0:(lx - 1),pars1[1:3],ddep)
       #lavec = lambdamu[[1]]
       #muvec = lambdamu[[2]]
       #m1 = lavec[1:lx,2:(lx+1)] * nx1[1:lx,2:(lx+1)]
       #m2 = muvec[3:(lx+2),2:(lx+1)] * nx1[3:(lx+2),2:(lx+1)]
       #ma = lavec[2:(lx+1),2:(lx+1)] + muvec[2:(lx+1),2:(lx+1)]
       #m3 = ma * nx1[2:(lx+1),2:(lx+1)]
       #m4 = lavec[2:(lx+1),1:lx] * nx2[2:(lx+1),1:lx]
       #m5 = muvec[2:(lx+1),3:(lx+2)] * nx2[2:(lx+1),3:(lx+2)]
       #m6 = ma * nx2[2:(lx+1),2:(lx+1)]
       probs = matrix(0,lx,lx)
       # probs[1,1] = probability of extinction of both lineages
       # sum(probs[1:lx,1]) = probability of extinction of second lineage
       probs[2,2] = 1 # clade M starts with two species
       # STEP 1: integrate from tcrown to tinn
       dim(probs) = c(lx*lx,1)
       #y = deSolve::ode(probs,c(tcrown,tinn),dd_logliknorm_rhs2,list(m1,m2,m3,m4,m5,m6),rtol = reltol,atol = abstol, method = "ode45")
       #probs = y[2,2:(lx * lx + 1)]
       #dim(probs) = c(lx,lx)
       #probs[1,1:lx] = 0                                                                           
       #probs[1:lx,1] = 0
       ly = lx^2
       #probs = rep(0,ly)
       #probs[lx + 2] = 1
       probs = dd_loglik_M2(pars = pars1[1:3],lx = lx,ddep = ddep,tt = abs(tinn - tcrown),p = probs)
       dim(probs) = c(lx,lx)
       probs[1,1:lx] = 0
       probs[1:lx,1] = 0
       #probs[1:lx] = 0                  
       #probs[seq(lx + 1,ly,by = lx)] = 0
       # STEP 2: transformation at tinn
       #dim(probs) = c(lx,lx)
       #nx1a = nx1[2:(lx+1),2:(lx+1)]
       #nx2a = nx2[2:(lx+1),2:(lx+1)]
       #probs = probs * nx1a/(nx1a+nx2a)
       #probs = rbind(probs[2:lx,1:lx], rep(0,lx))
       #dim(probs) = c(lx * lx,1)
       auxM1 = rep(0:(lx - 1),times = lx) + rep(0:(lx - 1),each = lx)
       probs = probs * rep(0:(lx - 1),times = lx)/auxM1
       dim(probs) = c(lx,lx)
       probs = rbind(probs[2:lx,1:lx], rep(0,lx))
       dim(probs) = c(lx * lx,1)
       #probs = c(probs[-c(seq(lx + 1,ly,by = lx))],rep(0,lx))
       # STEP 3: integrate from tinn to tpres
       #y = deSolve::ode(probs,c(tinn,tpres),dd_logliknorm_rhs2,list(m1,m2,m3,m4,m5,m6),rtol = reltol,atol = abstol, method = "ode45")
       #probs = y[2,2:(lx * lx + 1)]
       probs = dd_loglik_M2(pars = pars1[1:3],lx = lx,ddep = ddep,tt = abs(tpres - tinn),p = probs)
       dim(probs) = c(lx,lx)      
       PM12 = sum(probs[2:lx,2:lx])
       PM2 = sum(probs[1,2:lx])
       logliknorm = log(2) + (cond == 1) * log(PM12 + PS * PM2) + (cond == 4) * (log(PM12 + PM2) + log(PS))
    }
    if(length(m) > 1)
    {
       Sv = c(S1,S2)
    } else {
       Sv = S
    }
    loglik = loglik - logliknorm - sum(lgamma(Sv + m + 1) - lgamma(Sv + 1) - lgamma(m + 1))
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

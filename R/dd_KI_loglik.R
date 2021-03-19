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
#' @param pars1 Vector of parameters: \cr \cr
#' \code{pars1[1]} corresponds to
#' lambda_M (speciation rate) of the main clade \cr
#' \code{pars1[2]} corresponds
#' to mu_M (extinction rate) of the main clade \cr 
#' \code{pars1[3]} corresponds
#' to K_M (clade-level carrying capacity) of the main clade \cr
#' \code{pars1[4]} corresponds to lambda_S (speciation rate) of the subclade \cr
#' \code{pars1[5]} corresponds to mu_S (extinction rate) of the subclade \cr
#' \code{pars1[6]} corresponds to K_S (clade-level carrying capacity) of the
#' subclade \cr
#' \code{pars1[7]} corresponds to t_d (the time of decoupling)
#' @param pars2 Vector of model settings: \cr \cr
#' \code{pars2[1]} sets the
#' maximum number of species for which a probability must be computed.  This
#' must be larger than 1 + missnumspec + length(brts). \cr \cr
#' \code{pars2[2]} sets the model of diversity-dependence: \cr
#' - \code{pars2[2] == 1} linear
#' dependence in speciation rate with parameter K (= diversity where speciation
#' = extinction)\cr
#' - \code{pars2[2] == 1.3} linear dependence in speciation
#' rate with parameter K' (= diversity where speciation = 0)\cr
#' - \code{pars2[2] == 2} exponential dependence in speciation rate with
#' parameter K (= diversity where speciation = extinction)\cr 
#' - \code{pars2[2] == 2.1} variant of exponential dependence in speciation rate with offset at
#' infinity\cr
#' - \code{pars2[2] == 2.2} 1/n dependence in speciation rate\cr
#' - \code{pars2[2] == 2.3} exponential dependence in speciation rate with
#' parameter x (= exponent)\cr
#' - \code{pars2[2] == 3} linear dependence in
#' extinction rate \cr
#' - \code{pars2[2] == 4} exponential dependence in
#' extinction rate \cr
#' - \code{pars2[2] == 4.1} variant of exponential
#' dependence in extinction rate with offset at infinity\cr
#' - \code{pars2[2] == 4.2} 1/n dependence in extinction rate\cr\cr
#' \code{pars2[3]} sets the
#' conditioning: \cr
#' - \code{pars2[3] == 0} no conditioning (or just crown age) \cr
#' - \code{pars2[3] == 1} conditioning on non-extinction of the phylogeny \cr
#' - \code{pars2[3] == 2} conditioning on number of species and crown age;
#' not yet implemented \cr
#' - \code{pars2[3] == 3} conditioning on number of species only;
#' not yet implemented \cr
#' - \code{pars2[3] == 4} conditioning on survival of the subclade \cr
#' - \code{pars2[3] == 5} conditioning on survival of all subclades
#' and of both crown lineages in the main clade. This assumes that subclades
#' that have already shifted do not undergo another shift, i.e. shifts only
#' occur in the main clade. \cr \cr
#' \code{pars2[4]} Obsolete. \cr \cr
#' \code{pars2[5]} sets whether the parameters and likelihood should be
#' shown on screen (1) or not (0) \cr \cr
#' \code{pars2[6]} sets whether the first data point is stem age (1) or crown age (2) \cr\cr
#' \code{pars2[7]} sets whether the old (incorrect) likelihood should be used (0),
#' or whether the new corrected likelihood should be used (1).
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
dd_KI_loglik <- function(pars1,
                         pars2,
                         brtsM,
                         brtsS,
                         missnumspec,
                         methode = 'lsoda')
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
  #pars2 <- as.data.frame(rbind(pars2))
  names(pars2) = c('lx','ddep','cond','t_split','verbose','soc','corr')
  tinn = -abs(pars1[7])
  soc = pars2[6]
  #pars2 <- as.vector(pars2)
  # order branching times
  brts = -sort(abs(c(brtsM,brtsS)),decreasing = TRUE)
  if(sum(brts == 0) == 0)
  { 
    brts[length(brts) + 1] = 0
  }
  S = length(brts) + (soc - 2)
  brtsM = -sort(abs(brtsM),decreasing = TRUE)
  brtsS = -sort(abs(brtsS),decreasing = TRUE)
  tinn <- -abs(pars1[7])
  if(!any(brtsM == 0))
  { 
    brtsM = c(brtsM,0)
  }
  if(!any(brtsS == 0))
  { 
    brtsS = c(brtsS,0)
  }
  if(!any(brtsS == tinn))
  { 
    brtsS = c(tinn,brtsS)
  }
  
  # avoid coincidence of branching time and key innovation time
  if(sum(abs(brtsM - tinn) < 1E-14) == 1) { tinn = tinn - 1E-8 }
  
  ka = sum(brtsM < tinn)
  brtsMtinn <- sort(c(tinn,brtsM),decreasing = FALSE)
  kvec <- c(soc:(soc + ka - 1),
            (soc + ka - 2):(soc + length(brtsMtinn) - 4),
            soc + length(brtsMtinn) - 4)
  brts_kM <- rbind(brtsMtinn,kvec)
  
  kvec <- c(1:(length(brtsS)-1),length(brtsS)-1)
  brts_kS <- rbind(brtsS,kvec)
  
  brts_k_list <- list(brts_kM,brts_kS)
  pars1_list <- list(pars1[1:3],pars1[4:6])
  missnumspec_list <- create_missnumspec_list(missnumspec)
  are_pars_impossible <- check_for_impossible_pars_KI(S_list = create_S_list(brts_k_list,soc),
                                                      missnumspec_list = missnumspec_list,
                                                      pars1_list = pars1_list,
                                                      pars2 = pars2)
  if(are_pars_impossible || tinn >= 0 || tinn <= min(brtsM))
  {
    loglik <- -Inf
    return(loglik)
  }
  
  loglik <- dd_multiple_KI_loglik(pars1_list = pars1_list,
                                  pars2 = pars2,
                                  brts_k_list = brts_k_list,
                                  missnumspec_list = missnumspec_list,
                                  reltol = 1e-14,
                                  abstol = 1e-16,
                                  methode = methode)
  
  if(pars2[5])
  {
    s1 = sprintf('Parameters: %f %f %f %f %f %f %f, ',pars1[1],pars1[2],pars1[3],pars1[4],pars1[5],pars1[6],pars1[7])
    s2 = sprintf('Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    utils::flush.console()
  }
  loglik = as.numeric(loglik)
  if(is.nan(loglik) || is.na(loglik))
  {
    loglik = -Inf
  }
  return(loglik)
}   

shift <- function(probs,pars2,k1)
{
  lp <- length(probs)
  if(pars2[7] == 1)
  {
    probs <- probs * 1/(k1 + (0:(lp - 1)))
  } else
    if(pars2[7] == 1.5)
    {
      probs <- probs * k1/(k1 + (0:(lp - 1)))
    } else
      if(pars2[7] == 2)
      {
        probs <- probs * (0:(lp - 1))/(k1 + (0:(lp - 1))) 
      } else
        if(pars2[7] == 3)
        {
          probs[1:(lp - 1)] <- probs[2:lp] * (1:(lp - 1))/(k1 + (1:(lp - 1))) 
          probs[lp] <- 0
        }
  return(probs)
}

dd_KI_loglik_partial <- function(brts_k,
                                 m,
                                 cp,
                                 pars1,
                                 pars2,
                                 lx,
                                 reltol = 1e-14,
                                 abstol = 1e-16,
                                 methode,
                                 log = TRUE)
{  
  ddep <- pars2[2]
  verbose <- pars2[5]
  for(i in 2:length(brts_k[1,]))
  {
    loglik <- cp[[1]]
    probs <- cp[[2]]  
    t1 <- brts_k[1,i - 1]
    t2 <- brts_k[1,i]
    k1 <- brts_k[2,i - 1]
    k2 <- brts_k[2,i]
    y <- dd_integrate(probs,
                      c(t1,t2),
                      'dd_loglik_rhs',
                      c(pars1,k1,ddep),
                      rtol = reltol,
                      atol = abstol,
                      method = methode)
    probs <- y[2,2:(lx + 1)]
    if(k2 > k1)
    {
      probs <- flavec(ddep,
                      pars1[1],
                      pars1[2],
                      pars1[3],
                      r = 0,
                      lx,
                      k1) * probs # speciation event
    } else if(k2 <= k1 && t2 != 0)
    {
      probs <- shift(probs = probs,
                     pars2 = pars2,
                     k1 = k1)
    }
    cp <- check_probs(loglik = loglik,
                      probs = probs,
                      verbose = verbose)
  }
  loglik <- cp[[1]]
  probs <- cp[[2]]
  if(log == TRUE)
  {
    loglik <- loglik + log(probs[1 + m])
  } else 
  {
    loglik <- exp(loglik) * probs[1 + m] 
  }
  return(loglik)
}

convolve_logliks <- function(missnumspec_list,loglik_list)
{
  if(length(missnumspec_list[[1]]) == 1)
  {
    loglik <- sum(unlist(loglik_list))
  } else
  {
    loglik <- 0
    llist <- length(loglik_list)
    loglikdelta <- list()
    convolve_lik <- 1
    for(i in 1:llist)
    {
      loglikmax <- max(loglik_list[[i]])
      loglikdelta[[i]] <- loglik_list[[i]] - loglikmax
      loglik <- loglik + loglikmax
      convolve_lik <- conv(convolve_lik,exp(loglikdelta[[i]]))
    }
    loglik <- loglik + log(convolve_lik[max(missnumspec_list[[1]]) + 1])
  }
  if(is.nan(loglik) || is.na(loglik))
  {
    loglik <- -Inf
  }   
  return(loglik)
}

dd_KI_logliknorm <- function(brts_k_list,
                             pars1_list,
                             loglik,
                             cond,
                             ddep,
                             lx_list,
                             reltol = 1e-14,
                             abstol = 1e-16,
                             methode)
{
  if(cond == 0 || loglik == -Inf)
  {
    logliknorm = 0
  } else {
    if(length(brts_k_list) > 2 && cond == 1)
    {
      stop('Conditioning on survival not implemented for more than 1 shift')
    }
    if(cond == 2 || cond == 3)
    {
      stop('Conditioning on number of species not implemented')
    }
    tcrown <- min(brts_k_list[[1]][1,])
    ln <- length(brts_k_list[[1]][1,])
    kinn <- 1 + min(which(brts_k_list[[1]][2,2:ln] - brts_k_list[[1]][2,1:(ln - 1)] < 0)) 
    tinn <- brts_k_list[[1]][1,kinn]
    tpres <- 0
    # COMPUTE NORMALIZATION
    # compute survival probability of clade S
    lx = lx_list[[2]]
    nx = -1:lx
    lambdamu_nk <- lambdamu(nx,pars = c(pars1_list[[2]],0),ddep = ddep)
    lavec <- lambdamu_nk[[1]]
    muvec <- lambdamu_nk[[2]]
    probs = rep(0,lx) # probs[1] = extinction probability
    probs[2] = 1 # clade S starts with one species
    if(methode != 'analytical')
    {
      m1 = lavec[1:lx] * nx[1:lx]
      m2 = muvec[3:(lx + 2)] * nx[3:(lx + 2)]
      m3 = (lavec[2:(lx + 1)] + muvec[2:(lx + 1)]) * nx[2:(lx + 1)]

      if (startsWith(methode, 'odeint::')) {
        probs = .Call('dd_logliknorm1_odeint', probs, c(tinn,tpres), c(m1,m2,m3), abstol, reltol, methode)
      }
      else {
        y = deSolve::ode(probs,c(tinn,tpres),dd_logliknorm_rhs1,c(m1,m2,m3),rtol = reltol,atol = abstol,method = methode)
        probs = y[2,2:(lx+1)]
      }

    } else
    {
      probs = dd_loglik_M(pars1_list[[2]],lx,0,ddep,tt = abs(tpres - tinn),probs)
    }
    PS = 1 - probs[1]
    
    # compute survival probability of clade M
    lx = lx_list[[1]]
    n <- -1:lx
    nx1 = rep(n,lx + 2)
    dim(nx1) = c(lx + 2,lx + 2) # row index = number of species in first group 
    nx2 = t(nx1) # column index = number of species in second group
    nxt = nx1 + nx2
    lambdamu_n <- lambdamu2(nxt,pars1_list[[1]],ddep)
    lavec <- lambdamu_n[[1]]
    muvec <- lambdamu_n[[2]]
    probs = matrix(0,lx,lx)
    # probs[1,1] = probability of extinction of both lineages
    # sum(probs[1:lx,1]) = probability of extinction of second lineage
    probs[2,2] = 1 # clade M starts with two species
    # STEP 1: integrate from tcrown to tinn
    dim(probs) = c(lx,lx)
    if(methode != 'analytical')
    {
      m1 = lavec[1:lx,2:(lx+1)] * nx1[1:lx,2:(lx+1)]
      m2 = muvec[3:(lx+2),2:(lx+1)] * nx1[3:(lx+2),2:(lx+1)]
      ma = lavec[2:(lx+1),2:(lx+1)] + muvec[2:(lx+1),2:(lx+1)]
      m3 = ma * nx1[2:(lx+1),2:(lx+1)]
      m4 = lavec[2:(lx+1),1:lx] * nx2[2:(lx+1),1:lx]
      m5 = muvec[2:(lx+1),3:(lx+2)] * nx2[2:(lx+1),3:(lx+2)]
      m6 = ma * nx2[2:(lx+1),2:(lx+1)]
      if (startsWith(methode, "odeint::")) {
        probs = .Call('dd_logliknorm2_odeint', probs, c(tcrown,tinn), list(m1,m2,m3,m4,m5,m6), reltol, abstol, methode)
      }
      else {
        dim(probs) = c(lx*lx,1)
        y = deSolve::ode(probs,c(tcrown,tinn),dd_logliknorm_rhs2,list(m1,m2,m3,m4,m5,m6),rtol = reltol,atol = abstol, method = "ode45") 
        probs = y[2,2:(lx * lx + 1)]
      }
    } else
    {
      probs = dd_loglik_M2(pars = pars1_list[[1]],lx = lx,ddep = ddep,tt = abs(tinn - tcrown),p = probs)
    }
    dim(probs) = c(lx,lx)
    probs[1,1:lx] = 0
    probs[1:lx,1] = 0
    # STEP 2: transformation at tinn
    nx1a = nx1[2:(lx+1),2:(lx+1)]
    nx2a = nx2[2:(lx+1),2:(lx+1)]
    probs = probs * nx1a/(nx1a+nx2a)
    probs = rbind(probs[2:lx,1:lx], rep(0,lx))
    # STEP 3: integrate from tinn to tpres
    if(methode != 'analytical')
    {
      if (startsWith(methode, "odeint::")) {
        probs = .Call('dd_logliknorm2_odeint', probs, c(tinn, tpres), list(m1,m2,m3,m4,m5,m6), reltol, abstol, 'odeint::runge_kutta_fehlberg78')
      }
      else {
        dim(probs) = c(lx*lx,1)
        y = deSolve::ode(probs,c(tinn, tpres),dd_logliknorm_rhs2,list(m1,m2,m3,m4,m5,m6),rtol = reltol,atol = abstol, method = "ode45") 
        probs = y[2,2:(lx * lx + 1)]
        dim(probs) = c(lx,lx)
      }
    } else
    {
      probs = dd_loglik_M2(pars = pars1_list[[1]],lx = lx,ddep = ddep,tt = abs(tpres - tinn),p = probs)
    }
    dim(probs) = c(lx,lx)
    PM12 = sum(probs[2:lx,2:lx])
    PM2 = sum(probs[1,2:lx])
    #print(log(2) + log(PM12 + PS * PM2))
    #print(log(2) + (log(PM12 + PM2) + log(PS)))
    #print(log(2) + (log(PM12) + log(PS)))
    logliknorm = log(2) + (cond == 1) * log(PM12 + PS * PM2) +
      (cond == 4) * (log(PM12 + PM2) + log(PS)) +
      (cond == 5) * (log(PM12) + log(PS))
  }
  return(logliknorm)
}

create_lx_list <- function(lmax,
                           ddep,
                           missnumspec_list,
                           pars1_list)
{ 
  lx_list <- list()
  for(i in 1:length(pars1_list))
  {
    if(ddep == 1)
    {
      lx_list[[i]] <- min(max(1 + max(missnumspec_list[[i]]),1 + ceiling(pars1_list[[i]][1]/(pars1_list[[i]][1] - pars1_list[[i]][2]) * pars1_list[[i]][3])),ceiling(lmax))
    } else {
      if(ddep == 1.3)
      {
        lx_list[[i]] <- min(max(1 + max(missnumspec_list[[i]]),1 + ceiling(pars1_list[[i]][3])),ceiling(lmax))
      } else {
        lx_list[[i]] <- round(lmax)
      }
    }
  }
  return(lx_list)
}

create_missnumspec_list <- function(m)
{
  if(length(m) == 1)
  {
    missnumspec_list <- list(0:m,0:m)
  } else
  {
    missnumspec_list <- as.list(m)
  }
  return(missnumspec_list)
}

create_S_list <- function(brts_k_list,soc)
{
  S_list <- list()
  for(i in 1:length(brts_k_list))
  {
    S_list[[i]] <- brts_k_list[[i]][2,length(brts_k_list[[i]][2,])]
    names(S_list[[i]]) <- NULL
  }
  return(S_list)
}

check_for_impossible_pars_KI <- function(S_list,
                                         missnumspec_list,
                                         pars1_list,
                                         pars2)
{
  ddep <- pars2[2]
  loglik <- 0
  sumK <- 0
  sumSm <- max(missnumspec_list[[1]])
  for(i in 1:length(S_list))
  {
    are_pars_impossible <- check_for_impossible_pars(pars1_list[[i]],ddep)
    if(are_pars_impossible)
    {
      return(are_pars_impossible)
    }
    m <- 0
    if(length(missnumspec_list[[i]]) == 1)
    {
      m <- missnumspec_list[[i]]
    }
    are_pars_impossible <- 
      ((ddep == 1 && ceiling(pars1_list[[i]][1]/(pars1_list[[i]][1] - pars1_list[[i]][2]) * pars1_list[[i]][3]) < S_list[[i]] + m) ||
         (ddep == 1.3 && ceiling(pars1_list[[i]][3]) < S_list[[i]] + m))
    if(are_pars_impossible)
    {
      return(are_pars_impossible)
    }
    sumK <- sumK + ceiling(pars1_list[[i]][3])
    sumSm <- sumSm + S_list[[i]]
  }
  are_pars_impossible <- (ddep == 1.3 && sumK < sumSm)
  return(are_pars_impossible)
}

check_for_impossible_pars <- function(pars1,
                                      ddep)
{
  are_pars_impossible <- min(pars1) < 0
  if(are_pars_impossible)
  {
    cat('One or more parameters are negative.')
  } else
  { 
    are_pars_impossible <-   
      (pars1[2] == 0 && (ddep == 2 || ddep == 2.1 || ddep == 2.2)) ||
      (pars1[1] == 0 && (ddep == 4 || ddep == 4.1)) ||
      (pars1[1] <= pars1[2] && (ddep == 1 || ddep == 2 || ddep == 2.1 || ddep == 2.2 || ddep == 3 || ddep == 4 || ddep == 4.1))
    if(are_pars_impossible)
    {
      cat('These parameters are incompatible with the chosen model')
    }
  }
  return(are_pars_impossible)
}

#' Loglikelihood for diversity-dependent diversification models with multiple
#' decoupling (rate shift) events
#' 
#' This function computes loglikelihood of a diversity-dependent
#' diversification model for a given set of branching times and parameter
#' values where the diversity-dependent dynamics of subclades decouple from
#' the dynamics of main clades, potentially accompanied by a
#' shift in parameters.
#' @inheritParams dd_KI_loglik
#' @param pars1_list list of paramater sets one for each rate regime (subclade).
#' The parameters are: lambda (speciation rate), mu (extinction rate), and K
#' (clade-level carrying capacity).
#' @param brts_k_list list of matrices, one for each rate regime (subclade). Each
#' matrix has in the first row the branching times including the shift/decoupling time
#' and the present time (0) in negative time (i.e. 10 mya = -10). In the second row 
#' it has the number of lineages, i.e. starting at 2 for a phylogeny with a crown 
#' and increasing by one at each branching time and decreasing by one at each
#' decoupling/shift time. The last element is the same as the second last.
#' @param missnumspec_list list containing the number of missing species for each clade.
#' If only a single number m of missing species is known for the entire phylogeny, then each
#' element of the list should be 0:m. One can also create this from m using the function
#' create_missnumspec_list
#' @param reltol relative tolerance in integration of the ODE system, default at 1e-14
#' @param abstol tolerance tolerance in integration of the ODE system, default at 1e-16
dd_multiple_KI_loglik <- function(pars1_list,
                                  pars2,
                                  brts_k_list,
                                  missnumspec_list,
                                  reltol = 1e-14,
                                  abstol = 1e-16,
                                  methode = 'lsoda')
{
  lx_list <- create_lx_list(lmax = pars2[1],
                            ddep = pars2[2],
                            missnumspec_list = missnumspec_list,
                            pars1_list = pars1_list)
  llist <- length(pars1_list)
  if(llist != length(brts_k_list) || llist != length(missnumspec_list))
  {
    stop('The length of the lists of parameters, branching times and missing species is not the same.')
  }
  if(llist < 2)
  {
    stop('The length of the lists of parameters must be 2 or larger.')
  }
  loglik_list <- as.list(rep(0,llist))
  for(i in 1:llist)
  {
    probs <- rep(0,lx_list[[i]])
    probs[1] <- 1
    loglik_list[[i]] <- dd_KI_loglik_partial(brts_k = brts_k_list[[i]],
                                             m = missnumspec_list[[i]],
                                             cp = list(0,probs),
                                             pars1 = pars1_list[[i]],
                                             pars2 = pars2,
                                             lx = lx_list[[i]],
                                             reltol = reltol,
                                             abstol = abstol,
                                             methode = methode)
  }
  loglik <- convolve_logliks(missnumspec_list,loglik_list)
  if(pars2[3] != 5)
  {  
    logliknorm <- dd_KI_logliknorm(brts_k_list = brts_k_list,
                                   pars1_list = pars1_list,
                                   loglik = loglik,
                                   cond = pars2[3],
                                   ddep = pars2[2],
                                   lx_list = lx_list,
                                   reltol = reltol,
                                   abstol = abstol,
                                   methode = methode)
  } else
  {
    pars2[7] <- 3
    logliknorm <- dd_multiple_KI_logliknorm(brts_k_list = brts_k_list,
                                            pars1_list = pars1_list,
                                            pars2 = pars2,
                                            loglik = loglik,
                                            lx_list = lx_list,
                                            reltol = 1e-14,
                                            abstol = 1e-16,
                                            methode = methode)
  }
  loglik <- loglik - logliknorm
  S_list <- create_S_list(brts_k_list = brts_k_list,
                          soc = pars2[6])
  Sv <- unlist(S_list)
  m <- unlist(missnumspec_list)
  if(length(missnumspec_list[[1]]) > 1)
  {
    m <- max(missnumspec_list[[1]])
    Sv <- sum(Sv)
  }
  loglik <- loglik - sum(lgamma(Sv + m + 1) - lgamma(Sv + 1) - lgamma(m + 1))
  return(loglik)
}

dd_multiple_KI_logliknorm <- function(brts_k_list,
                                      pars1_list,
                                      pars2,
                                      loglik,
                                      lx_list,
                                      reltol = 1e-14,
                                      abstol = 1e-16,
                                      methode)
{
  if(pars2[3] == 0 || loglik == -Inf)
  {
    logliknorm <- 0
  } else
  {
    llist <- length(pars1_list)
    loglik_list <- as.list(rep(0,llist))
    logliknorm <- 0
    for(i in 1:llist)
    {
      probs <- rep(0,lx_list[[i]])
      probs[1] <- 1
      tstart <- min(brts_k_list[[i]][1,])
      ln <- length(brts_k_list[[i]][1,])
      kinn <- 1 + (which(brts_k_list[[i]][2,2:ln] - brts_k_list[[i]][2,1:(ln - 1)] < 0)) 
      tinn <- brts_k_list[[i]][1,kinn]
      tpres <- 0
      brts_k <- rbind(c(tstart,tinn,tpres),rep(brts_k_list[[i]][2,1],length(tinn) + 2))
      lp <- length(probs)
      if(brts_k_list[[i]][2,1] == 1)
      {
        aux <- 1:lp
      }
      if(brts_k_list[[i]][2,1] == 2)
      {
        aux <- (2:(lp + 1)) * (3:(lp + 2))/6
      }
      loglik_list[[i]] <- dd_KI_loglik_partial(brts_k = brts_k,
                                               m = 0:(lp - 1),
                                               cp = list(0,probs),
                                               pars1 = pars1_list[[i]],
                                               pars2 = pars2,
                                               lx = lx_list[[i]],
                                               reltol = reltol,
                                               abstol = abstol,
                                               methode = methode,
                                               log = FALSE)/aux
      logliknorm <- logliknorm + log(sum(loglik_list[[i]]))
    }
  }
  return(logliknorm)
}
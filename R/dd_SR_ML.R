#' Maximization of the loglikelihood under a diversity-dependent
#' diversification model with a shift in the parameters
#' 
#' This function computes the maximum likelihood estimates of the parameters of
#' a diversity-dependent diversification model with shifting parameters at time
#' t = tshift for a given set of phylogenetic branching times.  It also outputs
#' the corresponding loglikelihood that can be used in model comparisons.
#' 
#' The output is a dataframe containing estimated parameters and maximum
#' loglikelihood. The computed loglikelihood contains the factor q! m!/(q + m)!
#' where q is the number of species in the phylogeny and m is the number of
#' missing species, as explained in the supplementary material to Etienne et
#' al. 2012.
#' 
#' @param brts A set of branching times of a phylogeny, all positive
#' @param initparsopt The initial values of the parameters that must be
#' optimized
#' @param parsfix The values of the parameters that should not be optimized
#' @param idparsopt The ids of the parameters that must be optimized, e.g. 1:7
#' for all parameters.  The ids are defined as follows: \cr id == 1 corresponds
#' to lambda (speciation rate) before the shift \cr id == 2 corresponds to mu
#' (extinction rate) before the shift \cr id == 3 corresponds to K (clade-level
#' carrying capacity) before the shift \cr id == 4 corresponds to lambda
#' (speciation rate) after the shift \cr id == 5 corresponds to mu (extinction
#' rate) after the shift \cr id == 6 corresponds to K (clade-level carrying
#' capacity) after the shift \cr id == 7 corresponds to tshift (the time of
#' shift)
#' @param idparsfix The ids of the parameters that should not be optimized,
#' e.g. c(1,3,4,6) if lambda and K should not be optimized, but only mu. In
#' that case idparsopt must be c(2,5,7). The default is to fix all parameters
#' not specified in idparsopt.
#' @param idparsnoshift The ids of the parameters that should not shift; This
#' can only apply to ids 4, 5 and 6, e.g. idparsnoshift = c(4,5) means that
#' lambda and mu have the same values before and after tshift
#' @param res sets the maximum number of species for which a probability must
#' be computed, must be larger than 1 + length(brts)
#' @param ddmodel sets the model of diversity-dependence: \cr
#' ddmodel == 1 : linear dependence in speciation rate \cr
#' ddmodel == 2 : exponential dependence in speciation rate \cr
#' ddmodel == 2.1 : variant of exponential dependence in speciation rate with offset at infinity\cr ddmodel == 2.2 :
#' 1/n dependence in speciation rate\cr
#' ddmodel == 3 : linear dependence in extinction rate \cr
#' ddmodel == 4 : exponential dependence in extinction rate\cr
#' ddmodel == 4.1 : variant of exponential dependence in extinction rate
#' with offset at infinity\cr
#' ddmodel == 4.2 : 1/n dependence in extinction rate
#' with offset at infinity
#' @param missnumspec The number of species that are in the clade but missing
#' in the phylogeny
#' @param cond Conditioning: \cr
#' cond == 0 : no conditioning \cr
#' cond == 1 : conditioning on non-extinction of the phylogeny \cr
#' cond == 2 : conditioning on non-extinction of the phylogeny and on the total
#' number of extant taxa (including missing species) \cr
#' cond == 3 : conditioning on the total number of extant taxa (including missing
#' species) \cr
#' \cr Note: cond == 3 assumes a uniform prior on stem age, as is the standard
#' in constant-rate birth-death models, see e.g. D. Aldous & L. Popovic 2004.
#' Adv. Appl. Prob. 37: 1094-1115 and T. Stadler 2009. J. Theor. Biol. 261:
#' 58-66.
#' @param btorph Sets whether the likelihood is for the branching times (0) or
#' the phylogeny (1)
#' @param soc Sets whether stem or crown age should be used (1 or 2)
#' @param allbp Sets whether a search should be done with various initial
#' conditions, with tshift at each of the branching points (TRUE/FALSE)
#' @param tol Sets the tolerances in the optimization. Consists of: \cr reltolx
#' = relative tolerance of parameter values in optimization \cr reltolf =
#' relative tolerance of function value in optimization \cr abstolx = absolute
#' tolerance of parameter values in optimization
#' @param maxiter Sets the maximum number of iterations in the optimization
#' @param changeloglikifnoconv if TRUE the loglik will be set to -Inf if ML
#' does not converge
#' @param optimmethod Method used in optimization of the likelihood. Current
#' default is 'subplex'. Alternative is 'simplex' (default of previous
#' versions)
#' @param num_cycles the number of cycles of opimization. If set at Inf, it will
#' do as many cycles as needed to meet the tolerance set for the target function.
#' @param methode The method used to solve the master equation, default is
#' 'analytical' which uses matrix exponentiation; alternatively numerical ODE
#' solvers can be used, such as 'lsoda' or 'ode45'. These were used in the
#' package before version 3.1.
#' @param verbose Show the parameters and loglikelihood for every call to the
#' loglik function 
#' @return \item{lambda_1}{ gives the maximum likelihood estimate of lambda
#' before the shift} \item{mu_1}{ gives the maximum likelihood estimate of mu
#' before the shift} \item{K_1}{ gives the maximum likelihood estimate of K
#' before the shift} \item{lambda_2}{ gives the maximum likelihood estimate of
#' lambda after the shift} \item{mu_2}{ gives the maximum likelihood estimate
#' of mu after the shift} \item{K_2}{ gives the maximum likelihood estimate of
#' K after the shift} \item{t_shift}{ gives the time of the shift}
#' \item{loglik}{ gives the maximum loglikelihood} \item{df}{ gives the number
#' of estimated parameters, i.e. degrees of feedom} \item{conv}{ gives a
#' message on convergence of optimization; conv = 0 means convergence}
#' @note The optimization may get trapped in local optima. Try different
#' starting values to search for the global optimum.
#' @author Rampal S. Etienne & Bart Haegeman
#' @seealso \code{\link{dd_SR_loglik}}, \code{\link{dd_ML}},
#' \code{\link{dd_KI_ML}},
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#' 
#' cat("This will estimate parameters for a sets of branching times brts without conditioning.\n")
#' cat("The tolerance of the optimization is set ridiculously high to make runtime fast.\n")
#' cat("In real applications, use the default or more stringent settings for tol.\n")
#' brts = 1:10
#' dd_SR_ML(brts = brts, initparsopt = c(0.4581, 1E-6, 17.69, 11.09, 8.9999), idparsopt = c(1:3,6,7),
#'          idparsfix = NULL, parsfix = NULL, idparsnoshift = c(4,5), cond = 0,
#'          tol = c(1E-1,1E-1,1E-1),optimmethod = 'simplex'
#' )
#' 
#' @export dd_SR_ML
dd_SR_ML = function(brts,
    initparsopt = c(0.5,0.1,2*(1+length(brts)+missnumspec),2*(1+length(brts)+missnumspec),max(brts)/2),
    parsfix = NULL,
    idparsopt = c(1:3,6:7),
    idparsfix = NULL,
    idparsnoshift = (1:7)[c(-idparsopt,(-1)^(length(idparsfix) != 0) * idparsfix)],
    res = 10*(1 + length(brts) + missnumspec),
    ddmodel = 1,
    missnumspec = 0,
    cond = 1,
    btorph = 1,
    soc = 2,
    allbp = FALSE,
    tol = c(1E-3, 1E-4, 1E-6),
    maxiter = 1000 * round((1.25)^length(idparsopt)),
    changeloglikifnoconv = FALSE,
    optimmethod = 'subplex',
    num_cycles = 1,
    methode = 'analytical',
    verbose = FALSE)
{
# brts = branching times (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# - idparsopt contains the ids of the parameters to be optimized, e.g. to optimize la, mu, K, K2 and tshift idparsopt = c(1,2,3,6,7)
# - initparsopt contains the starting values of the parameters to be optimized
# - idparsfix contains the ids of the parameters that are fixed and must not be optimized
# - parsfix contains the values of the fixed parameters
# - idparsnoshift contains the ids of the parameters la2, mu2 and K2 that do not shift, i.e. that need to be set equal to la, mu and K
# - pars[1] = la = (initial) speciation rate before shift
# - pars[2] = mu = extinction rate before shift
# - pars[3] = K = carrying capacity before shift
# - pars[4] = la2 = (initial) speciation rate after shift
# - pars[5] = mu2 = extinction rate after shift
# - pars[6] = K2 = carrying capacity after shift
# - pars[7] = tshift = time of shift
# - res = resolution of the method; res should be larger than the total number of species
# - ddmodel = diversity-dependent model,mode of diversity-dependence
#  . ddmodel == 1 : linear dependence in speciation rate with parameter K
#  . ddmodel == 1.3: linear dependence in speciation rate with parameter K'
#  . ddmodel == 2 : exponential dependence in speciation rate
#  . ddmodel == 2.1: variant with offset at infinity
#  . ddmodel == 2.2: 1/n dependence in speciation rate
#  . ddmodel == 2.3: exponential dependence in speciation rate with parameter x
#  . ddmodel == 3 : linear dependence in extinction rate
#  . ddmodel == 4 : exponential dependence in extinction rate
#  . ddmodel == 4.1: variant with offset at infinity
#  . ddmodel == 4.2: 1/n dependence in speciation rate
# - missnumspec = number of missing species    
# - cond = conditioning
#  . cond == 0 : no conditioning
#  . cond == 1 : conditioning on non-extinction of the phylogeny
#  . cond == 2 : conditioning on non-extinction of the phylogeny and on the total number of extant taxa (including missing species)
# - btorph = likelihood of branching times (0) or phylogeny (1), differ by a factor (S - 1)! where S is the number of extant species
# - allbp = optimize likelihood for fixed tshift at all bp (TRUE) or by letting tshift vary freely (FALSE)
# - tol = tolerance in optimization
#  . reltolx = relative tolerance of parameter values in optimization
#  . reltolf = relative tolerance of function value in optimization
#  . abstolx = absolute tolerance of parameter values in optimization
# - maxiter = the maximum number of iterations in the optimization
# - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
# - optimmethod = 'subplex' (current default) or 'simplex' (default of previous versions)
# - methode = the method used in the numerical solving of the set of the ode's
  brts = sort(abs(as.numeric(brts)),decreasing = TRUE)
  #options(warn = -1)
  if(is.numeric(brts) == FALSE)
  {
    cat("The branching times should be numeric.\n")
    out2 = data.frame(row.names = "results",lambda_1 = -1, mu_1 = -1, K_1 = -1, lambda_2 = -1, mu_2 = -1, K_2 = -1, t_shift = -1, loglik = -1, df = -1, conv = -1)
  } else {
    idparsnoshift = sort(idparsnoshift)
    idpars = sort(c(idparsopt,idparsfix,idparsnoshift))
    if((prod(idpars == (1:7)) != 1) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)))
    {
      cat("The parameters to be optimized, fixed and not shifted are incoherent.\n")
      out2 = data.frame(row.names = "results",lambda_1 = -1, mu_1 = -1, K_1 = -1, lambda_2 = -1, mu_2 = -1, K_2 = -1, t_shift = -1, loglik = -1, df = -1, conv = -1)
    } else {
      namepars = c("la","mu","K","la2","mu2","K2","tshift")
      if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
      cat("You are optimizing",optstr,"\n")
      if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
      cat("You are fixing",fixstr,"\n")
      if(length(namepars[idparsnoshift]) == 0) { noshiftstr = "anything" } else { noshiftstr = namepars[idparsnoshift] }
      cat("You are not shifting",noshiftstr,"\n")
      cat("Optimizing the likelihood - this may take a while.","\n")
      utils::flush.console()
      trparsopt = initparsopt/(1 + initparsopt)
      trparsopt[which(initparsopt == Inf)] = 1
      trparsfix = parsfix/(1 + parsfix)
      trparsfix[which(parsfix == Inf)] = 1
      pars2 = c(res,ddmodel,cond,btorph,verbose,soc,tol,maxiter)
      optimpars = c(tol,maxiter)
      initloglik = dd_SR_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsnoshift = idparsnoshift,pars2 = pars2,brts = brts,missnumspec = missnumspec, methode = methode)
      cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
      if(initloglik == -Inf)
      {
        cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
        out2 = data.frame(row.names = "results",lambda_1 = -1, mu_1 = -1, K_1 = -1, lambda_2 = -1, mu_2 = -1, K_2 = -1, t_shift = -1, loglik = -1, df = -1, conv = -1)
      } else {
        #code up to DDD v1.6: out = optimx2(trparsopt,dd_SR_loglik_choosepar,hess=NULL,method = "Nelder-Mead",hessian = FALSE,control = list(maximize = TRUE,abstol = pars2[8],reltol = pars2[7],trace = 0,starttests = FALSE,kkt = FALSE),trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsnoshift = idparsnoshift,brts = brts,pars2 = pars2,missnumspec = missnumspec)
        out = optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = dd_SR_loglik_choosepar,trparsopt = trparsopt,idparsopt = idparsopt,trparsfix = trparsfix,idparsfix = idparsfix,idparsnoshift = idparsnoshift,pars2 = pars2,brts = brts,missnumspec = missnumspec, methode = methode, num_cycles = num_cycles)
        if(out$conv != 0)
        {
          cat("Optimization has not converged. Try again with different initial values.\n")
          out2 = data.frame(row.names = "results",lambda_1 = -1, mu_1 = -1, K_1 = -1, lambda_2 = -1, mu_2 = -1, K_2 = -1, t_shift = -1, loglik = -1, df = -1, conv = unlist(out$conv))
        } else {
          MLtrpars = as.numeric(unlist(out$par))
          MLpars = MLtrpars/(1 - MLtrpars)
          out$par = list(MLpars)
          MLpars1 = rep(0,7)
          MLpars1[idparsopt] = MLpars
          ML = as.numeric(unlist(out$fvalues))
          if(sum(idparsfix == 7) == 0 && allbp == TRUE)
          { 
            idparsopt1 = idparsopt[1:(length(idparsopt) - 1)]
            idparsfix1 = c(idparsfix,7)
            for(bp in 2:length(brts))
            {
              for(ba in seq(-1,1,2))
              {
                initparsopt1 = initparsopt[1:length(idparsopt1)]
                parsfix1 = c(idparsfix,brts[bp] + ba * 1E-8)
                trparsopt1 = initparsopt1/(1 + initparsopt1)
                trparsfix1 = parsfix1/(1 + parsfix1)
                #code up to DDD v1.6: out = optimx2(trparsopt1,dd_SR_loglik_choosepar,hess=NULL,method = "Nelder-Mead",hessian = FALSE,control = list(maximize = TRUE,abstol = 1E-10,reltol = pars2[2],trace = 0,starttests = FALSE,kkt = FALSE),trparsfix = trparsfix1,idparsopt = idparsopt1,idparsfix = idparsfix1,idparsnoshift = idparsnoshift,brts = brts,pars2 = pars2,missnumspec = missnumspec)
                out = optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = dd_SR_loglik_choosepar,trparsopt = trparsopt1,idparsopt = idparsopt1,trparsfix = trparsfix1,idparsfix = idparsfix1,idparsnoshift = idparsnoshift,pars2 = pars2,brts = brts,missnumspec = missnumspec, methode = methode, num_cycles = num_cycles)
                if(as.numeric(out$fvalues) > ML)
                {
                  MLtrpars = as.numeric(unlist(out$par))
                  MLpars = MLtrpars/(1-MLtrpars)
                  out$par = list(MLpars)
                  ML = as.numeric(unlist(out$fvalues))
                }
              }
            }
          }
          if(length(idparsfix) != 0) {MLpars1[idparsfix] = parsfix }
          if(length(idparsnoshift) != 0) { MLpars1[idparsnoshift] = MLpars1[idparsnoshift - 3] }
          if(MLpars1[3] > 10^7){ MLpars1[3] = Inf}
          if(MLpars1[6] > 10^7){ MLpars1[6] = Inf}
          s1 = sprintf('Maximum likelihood parameter estimates: %f %f %f %f %f %f %f',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7])
          s2 = sprintf('Maximum loglikelihood: %f',ML)
          cat(paste("\n",s1,"\n",s2,"\n\n",sep = ''))
          out$par = list(MLpars1)
          out$fvalues = list(ML)
          out2 = data.frame(row.names = "results",lambda_1 = MLpars1[1],mu_1 = MLpars1[2],K_1 = MLpars1[3],lambda_2 = MLpars1[4],mu_2 = MLpars1[5],K_2 = MLpars1[6],t_shift = MLpars1[7],loglik = ML,df = length(initparsopt),conv = unlist(out$conv))
          if(out2$conv != 0 & changeloglikifnoconv == T) { out2$loglik = -Inf }
        }
      }
    }
  }
  return(invisible(out2))
}

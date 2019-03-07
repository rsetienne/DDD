#' Maximization of the loglikelihood under the diversity-independent, possibly
#' time-dependent diversification model
#' 
#' This function computes the maximum likelihood estimates of the parameters of
#' a diversity-independent diversification model for a given set of
#' phylogenetic branching times.  It also outputs the corresponding
#' loglikelihood that can be used in model comparisons.
#' 
#' The output is a dataframe containing estimated parameters and maximum
#' loglikelihood. The computed loglikelihood contains the factor q! m! / (q +
#' m)! where q is the number of species in the phylogeny and m is the number of
#' missing species, as explained in the supplementary material to Etienne et
#' al. 2012.
#' 
#' @param brts A set of branching times of a phylogeny, all positive
#' @param initparsopt The initial values of the parameters that must be
#' optimized
#' @param idparsopt The ids of the parameters that must be optimized, e.g. 1:3
#' for intrinsic speciation rate, extinction rate and carrying capacity.  The
#' ids are defined as follows: \cr id == 1 corresponds to lambda0 (speciation
#' rate) \cr id == 2 corresponds to mu0 (extinction rate) \cr id == 3
#' corresponds to lamda1 (parameter controlling decline in speciation rate with
#' time) \cr id == 4 corresponds to mu1 (parameter controlling decline in
#' extinction rate with time)
#' @param idparsfix The ids of the parameters that should not be optimized,
#' e.g. c(1,3) if lambda0 and lambda1 should not be optimized, but only mu0 and
#' mu1. In that case idparsopt must be c(2,4). The default is to fix all
#' parameters not specified in idparsopt.
#' @param parsfix The values of the parameters that should not be optimized
#' @param missnumspec The number of species that are in the clade but missing
#' in the phylogeny
#' @param tdmodel Sets the model of time-dependence: \cr
#' tdmodel == 0 : constant speciation and extinction rates \cr
#' tdmodel == 1 : speciation and/or extinction rate is exponentially declining
#' with time \cr
#' tdmodel == 2 : stepwise decline in speciation rate as in diversity-dependence
#' without extinction \cr
#' tdmodel == 3 : decline in speciation rate following deterministic logistic
#' equation for ddmodel = 1 \cr
#' tdmodel == 4 : decline in speciation rate such that the expected number of species matches with
#' that of ddmodel = 1 with the same mu
#' @param cond Conditioning: \cr
#' cond == 0 : conditioning on stem or crown age \cr
#' cond == 1 : conditioning on stem or crown age and non-extinction of the
#' phylogeny \cr
#' cond == 2 : conditioning on stem or crown age and on the total
#' number of extant taxa (including missing species) \cr
#' cond == 3 : conditioning on the total number of extant taxa (including missing species)
#' @param btorph Sets whether the likelihood is for the branching times (0) or
#' the phylogeny (1)
#' @param soc Sets whether stem or crown age should be used (1 or 2)
#' @param tol Sets the tolerances in the optimization. Consists of: \cr
#' reltolx = relative tolerance of parameter values in optimization \cr
#' reltolf = relative tolerance of function value in optimization \cr
#' abstolx = absolute tolerance of parameter values in optimization
#' @param maxiter Sets the maximum number of iterations in the optimization
#' @param changeloglikifnoconv if TRUE the loglik will be set to -Inf if ML
#' does not converge
#' @param optimmethod Method used in optimization of the likelihood. Current
#' default is 'subplex'. Alternative is 'simplex' (default of previous
#' versions)
#' @param num_cycles the number of cycles of opimization. If set at Inf, it will
#' do as many cycles as needed to meet the tolerance set for the target function.
#' @param methode The method used to solve the master equation under tdmodel =
#' 4, default is 'lsoda'.
#' @return A dataframe with the following elements:\cr
#' \item{lambda0}{ gives the maximum likelihood estimate of lambda0}
#' \item{mu0}{ gives the maximum likelihood estimate of mu0}
#' \item{lambda1}{gives the maximum likelihood estimate of lambda1}
#' \item{mu1}{ gives the
#' maximum likelihood estimate of mu1}
#' \item{loglik}{ gives the maximum
#' loglikelihood}
#' \item{df}{ gives the number of estimated parameters, i.e.
#' degrees of feedom}
#' \item{conv}{ gives a message on convergence of optimization; conv = 0 means
#' convergence}
#' @author Rampal S. Etienne & Bart Haegeman
#' @seealso \code{\link{bd_loglik}}
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#' 
#' cat("Estimating parameters for a set of branching times brts with the default settings:")
#' brts = 1:20
#' bd_ML(brts = brts, cond = 1)
#' 
#' @export bd_ML
bd_ML = function(brts, initparsopt = c(0.1,0.05 * (tdmodel <= 1) + 10 * (length(brts) + missnumspec) * (tdmodel > 1)), idparsopt = c(1,2 + (tdmodel > 1)), idparsfix = (1:4)[-idparsopt], parsfix = rep(0,4)[idparsfix], missnumspec = 0, tdmodel = 0, cond = 1, btorph = 1, soc = 2, tol = c(1E-3, 1E-4, 1E-6), maxiter = 1000 * round((1.25)^length(idparsopt)), changeloglikifnoconv = FALSE, optimmethod = 'subplex',num_cycles = 1, methode = 'lsoda')
{
  options(warn = -1)
  brts = sort(abs(as.numeric(brts)),decreasing = TRUE)
  out2 = invisible(data.frame(lambda0 = -1,mu0 = -1,lambda1 = -1, mu1 = -1, loglik = -1, df = -1, conv = -1))
  if(is.numeric(brts) == FALSE)
  {
     cat("The branching times should be numeric.\n")
     return(out2)
  }
  idpars = sort(c(idparsopt,idparsfix))
  if((prod(idpars == (1:4)) != 1) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)))
  {
     cat("The parameters to be optimized and/or fixed are incoherent.\n")
     return(out2)
  }
  if(tdmodel == 0 & length(initparsopt) > 2)
  {
     cat("tdmodel = 0 only accepts a constant speciation rate and extinction rate.\n")
     return(out2)
  }
  namepars1 = c("lambda0","mu0","lambda1","mu1")
  namepars2 = c("lambda0","mu0","K","")
  if(tdmodel == 2 | tdmodel == 3) { namepars = namepars2 } else  { namepars = namepars1 }
  if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
  cat("You are optimizing",optstr,"\n")
  if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
  cat("You are fixing",fixstr,"\n")
  cat("Optimizing the likelihood - this may take a while.","\n")
  flush.console()
  trparsopt = initparsopt/(1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] = 1
  trparsfix = parsfix/(1 + parsfix)
  trparsfix[which(parsfix == Inf)] = 1
  pars2 = c(tdmodel,cond,btorph,0,soc,1000,tol,maxiter)
  optimpars = c(tol,maxiter)
  initloglik = bd_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,pars2 = pars2,brts = brts,missnumspec = missnumspec, methode = methode)
  cat("The loglikelihood for the inital parameter values is",initloglik,"\n")
  if(initloglik == -Inf)
  {
     cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
     return(out2)
  }
  out = optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = bd_loglik_choosepar,trparsopt = trparsopt,idparsopt = idparsopt,trparsfix = trparsfix,idparsfix = idparsfix,pars2 = pars2,brts = brts,missnumspec = missnumspec,methode = methode, num_cycles = num_cycles)
  if(out$conv != 0)
  {
     cat("Optimization has not converged. Try again with different initial values.\n")
     return(out2)
  }
  MLtrpars = as.numeric(unlist(out$par))
  MLpars = MLtrpars/(1-MLtrpars)
  MLpars1 = rep(0,4)
  MLpars1[idparsopt] = MLpars
  if(length(idparsfix) != 0) { MLpars1[idparsfix] = parsfix }
  ML = as.numeric(unlist(out$fvalues))
  out2 = data.frame(lambda0 = MLpars1[1], mu0 = MLpars1[2], lambda1 = MLpars1[3], mu1 = MLpars1[4], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
  if(out2$conv != 0 & changeloglikifnoconv == T) { out2$loglik = -Inf }
  s1 = sprintf('Maximum likelihood parameter estimates: lambda0: %f, mu0: %f, lambda1: %f, mu1: %f: ',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4])
  s2 = sprintf('Maximum loglikelihood: %f',ML)
  cat("\n",s1,"\n",s2,"\n")
  out2 = invisible(out2)
  return(out2)
}

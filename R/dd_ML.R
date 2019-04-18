initparsoptdefault = function(ddmodel,brts,missnumspec)
{
  if(ddmodel < 5)
  {
    return(c(0.2,0.1,2 * (length(brts) + missnumspec)^(ddmodel != 2.3)))
  } else {
    return(c(0.2,0.1,2 * (length(brts) + missnumspec),0.01))
  }
}

parsfixdefault = function(ddmodel,brts,missnumspec,idparsopt)
{
  if(ddmodel < 5)
  {
    return(c(0.2,0.1,2*(length(brts) + missnumspec))[-idparsopt])
  } else {
    return(c(0.2,0.1,2*(length(brts) + missnumspec),0)[-idparsopt])
  }
}



#' Maximization of the loglikelihood under a diversity-dependent
#' diversification model
#' 
#' This function computes the maximum likelihood estimates of the parameters of
#' a diversity-dependent diversification model for a given set of phylogenetic
#' branching times.  It also outputs the corresponding loglikelihood that can
#' be used in model comparisons.
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
#' ids are defined as follows: \cr id == 1 corresponds to lambda (speciation
#' rate) \cr id == 2 corresponds to mu (extinction rate) \cr id == 3
#' corresponds to K (clade-level carrying capacity) \cr id == 4 corresponds to
#' r (r = b/a where mu = mu_0 + b * N and lambda = lambda_0 - a * N) (This is
#' only available when ddmodel = 5)
#' @param idparsfix The ids of the parameters that should not be optimized,
#' e.g. c(1,3) if lambda and K should not be optimized, but only mu. In that
#' case idparsopt must be 2. The default is to fix all parameters not specified
#' in idparsopt.
#' @param parsfix The values of the parameters that should not be optimized
#' @param res Sets the maximum number of species for which a probability must
#' be computed, must be larger than 1 + length(brts)
#' @param ddmodel Sets the model of diversity-dependence: \cr
#' \code{ddmodel == 1} : linear dependence in speciation rate with parameter K (= diversity
#' where speciation = extinction)\cr
#' \code{ddmodel == 1.3} : linear dependence
#' in speciation rate with parameter K' (= diversity where speciation = 0)\cr
#' \code{ddmodel == 1.4} : positive diversity-dependence in speciation rate
#' with parameter K' (= diversity where speciation rate reaches half its
#' maximum); lambda = lambda0 * S/(S + K') where S is species richness\cr
#' \code{ddmodel == 1.5} : positive and negative dependence in speciation rate
#' with parameter K' (= diversity where speciation = 0); lambda = lambda0 *
#' S/K' * (1 - S/K') where S is species richness\cr
#' \code{ddmodel == 2} : exponential dependence in speciation rate with parameter
#' K (= diversity where speciation = extinction)\cr
#' \code{ddmodel == 2.1} : variant of exponential dependence in speciation rate
#' with offset at infinity\cr
#' \code{ddmodel == 2.2} : 1/n dependence in speciation rate\cr 
#' \code{ddmodel == 2.3} : exponential dependence in speciation rate with parameter x (=
#' exponent)\cr
#' \code{ddmodel == 3} : linear dependence in extinction rate \cr
#' \code{ddmodel == 4} : exponential dependence in extinction rate \cr
#' \code{ddmodel == 4.1} : variant of exponential dependence in extinction rate
#' with offset at infinity \cr
#' \code{ddmodel == 4.2} : 1/n dependence in extinction rate with offset at infinity \cr \code{ddmodel == 5} : linear
#' dependence in speciation and extinction rate \cr
#' @param missnumspec The number of species that are in the clade but missing
#' in the phylogeny
#' @param cond Conditioning: \cr
#' cond == 0 : conditioning on stem or crown age\cr
#' cond == 1 : conditioning on stem or crown age and non-extinction of the
#' phylogeny \cr
#' cond == 2 : conditioning on stem or crown age and on the total
#' number of extant taxa (including missing species) \cr
#' cond == 3 : conditioning on the total number of extant taxa (including missing species)\cr
#' Note: cond == 3 assumes a uniform prior on stem age, as is the standard
#' in constant-rate birth-death models, see e.g. D. Aldous & L. Popovic 2004.
#' Adv. Appl. Prob. 37: 1094-1115 and T. Stadler 2009. J. Theor. Biol. 261:
#' 58-66.
#' @param btorph Sets whether the likelihood is for the branching times (0) or
#' the phylogeny (1)
#' @param soc Sets whether stem or crown age should be used (1 or 2)
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
#' @return \item{lambda}{ gives the maximum likelihood estimate of lambda}
#' \item{mu}{ gives the maximum likelihood estimate of mu} \item{K}{ gives the
#' maximum likelihood estimate of K} \item{r}{ (only if ddmodel == 5) gives the
#' ratio of linear dependencies in speciation and extinction rates}
#' \item{loglik}{ gives the maximum loglikelihood} \item{df}{ gives the number
#' of estimated parameters, i.e. degrees of feedom} \item{conv}{ gives a
#' message on convergence of optimization; conv = 0 means convergence}
#' @author Rampal S. Etienne & Bart Haegeman
#' @seealso \code{\link{dd_loglik}}, \code{\link{dd_SR_ML}},
#' \code{\link{dd_KI_ML}},
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#' 
#' cat("Estimating the intrinsic speciation rate lambda and the carrying capacity K")
#' cat("for a fixed extinction rate of 0.1, conditioning on clade survival and two missing species:")
#' brts = 1:5
#' dd_ML(brts = brts,initparsopt = c(1.3078,7.4188), idparsopt = c(1,3), parsfix = 0.1,
#'       cond = 1, missnumspec = 2, tol = c(1E-3,1E-3,1E-4), optimmethod = 'simplex')
#' 
#' @export dd_ML
dd_ML = function(
  brts,
  initparsopt = initparsoptdefault(ddmodel,brts,missnumspec),
  idparsopt = 1:length(initparsopt),
  idparsfix = (1:(3 + (ddmodel == 5)))[-idparsopt],
  parsfix = parsfixdefault(ddmodel,brts,missnumspec,idparsopt),
  res = 10*(1+length(brts)+missnumspec),
  ddmodel = 1,
  missnumspec = 0,
  cond = 1,
  btorph = 1,
  soc = 2,
  tol = c(1E-3, 1E-4, 1E-6),
  maxiter = 1000 * round((1.25)^length(idparsopt)),
  changeloglikifnoconv = FALSE,
  optimmethod = 'subplex',
  num_cycles = 1,
  methode = 'analytical',
  verbose = FALSE)
{
  #options(warn = -1)
  brts = sort(abs(as.numeric(brts)),decreasing = TRUE)
  if(is.numeric(brts) == FALSE)
  {
    cat("The branching times should be numeric.\n")
    out2 = data.frame(lambda = -1,mu = -1,K = -1, loglik = -1, df = -1, conv = -1)
    if(ddmodel == 5) {out2 = data.frame(lambda = -1,mu = -1,K = -1, r = -1, loglik = -1, df = -1, conv = -1)}
  } else {
    idpars = sort(c(idparsopt,idparsfix))
    if((prod(idpars == (1:(3 + (ddmodel == 5)))) != 1) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)))
    {
      cat("The parameters to be optimized and/or fixed are incoherent.\n")
      out2 = data.frame(lambda = -1,mu = -1,K = -1, loglik = -1, df = -1, conv = -1)
      if(ddmodel == 5) {out2 = data.frame(lambda = -1,mu = -1,K = -1, r = -1, loglik = -1, df = -1, conv = -1)}
    } else {
      namepars = c("lambda","mu","K")
      if(ddmodel == 5) {namepars = namepars = c("lambda","mu","K","r")}
      if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
      cat("You are optimizing",optstr,"\n")
      if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
      cat("You are fixing",fixstr,"\n")
      cat("Optimizing the likelihood - this may take a while.","\n")
      utils::flush.console()
      trparsopt = initparsopt/(1 + initparsopt)
      trparsopt[which(initparsopt == Inf)] = 1
      trparsfix = parsfix/(1 + parsfix)
      trparsfix[which(parsfix == Inf)] = 1
      pars2 = c(res,ddmodel,cond,btorph,verbose,soc,tol,maxiter)
      optimpars = c(tol,maxiter)
      initloglik = dd_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,pars2 = pars2,brts = brts,missnumspec = missnumspec, methode = methode)
      cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
      utils::flush.console()
      if(initloglik == -Inf)
      {
        cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
        out2 = data.frame(lambda = -1,mu = -1,K = -1, loglik = -1, df = -1, conv = -1)
        if(ddmodel == 5) {out2 = data.frame(lambda = -1,mu = -1,K = -1, r = -1, loglik = -1, df = -1, conv = -1)}
      } else {
        #code up to DDD v1.6: out = optimx2(trparsopt,dd_loglik_choosepar,hess=NULL,method = "Nelder-Mead",hessian = FALSE,control = list(maximize = TRUE,abstol = pars2[8],reltol = pars2[7],trace = 0,starttests = FALSE,kkt = FALSE),trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,brts = brts, pars2 = pars2,missnumspec = missnumspec)
        #out = dd_simplex(trparsopt,idparsopt,trparsfix,idparsfix,pars2,brts,missnumspec)
        out = optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = dd_loglik_choosepar,trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,pars2 = pars2,brts = brts, missnumspec = missnumspec, methode = methode, num_cycles = num_cycles)
        if(out$conv != 0)
        {
          cat("Optimization has not converged. Try again with different initial values.\n")
          out2 = data.frame(lambda = -1,mu = -1,K = -1, loglik = -1, df = -1, conv = unlist(out$conv))
          if(ddmodel == 5) {out2 = data.frame(lambda = -1,mu = -1,K = -1, r = -1, loglik = -1, df = -1, conv = unlist(out$conv))}
        } else {
          MLtrpars = as.numeric(unlist(out$par))
          MLpars = MLtrpars/(1-MLtrpars)
          MLpars1 = rep(0,3)
          if(ddmodel == 5) {MLpars1 = rep(0,4)}
          MLpars1[idparsopt] = MLpars
          if(length(idparsfix) != 0) { MLpars1[idparsfix] = parsfix }
          if(MLpars1[3] > 10^7){MLpars1[3] = Inf}
          ML = as.numeric(unlist(out$fvalues))
          out2 = data.frame(lambda = MLpars1[1],mu = MLpars1[2],K = MLpars1[3], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
          s1 = sprintf('Maximum likelihood parameter estimates: lambda: %f, mu: %f, K: %f',MLpars1[1],MLpars1[2],MLpars1[3])
          if(ddmodel == 5)
          {
            s1 = sprintf('%s, r: %f',s1,MLpars1[4])
            out2 = data.frame(lambda = MLpars1[1],mu = MLpars1[2],K = MLpars1[3], r = MLpars1[4], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))   
          }
          if(out2$conv != 0 & changeloglikifnoconv == T) { out2$loglik = -Inf }
          s2 = sprintf('Maximum loglikelihood: %f',ML)
          cat(paste("\n",s1,"\n",s2,"\n",sep = ''))
        }
      }
    }
  }
   return(invisible(out2))
}

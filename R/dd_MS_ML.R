#' Maximization of the loglikelihood under a diversity-dependent
#' diversification model with decoupling of a subclade's diversication dynamics
#' from the main clade's dynamics
#' 
#' This function computes the maximum likelihood estimates of the parameters of
#' a diversity-dependent diversification model where the diversity-dependent
#' dynamics of an innovative subclade have different parameters from the
#' dynamics of the main clade from time t_d, but both are governed by the same
#' carrying capacity and experience each other's diversity.  Required isa given
#' set of phylogenetic branching times of main clade and subclade and the time
#' of splitting of the lineage that will form the subclade.  The function also
#' outputs the corresponding loglikelihood that can be used in model
#' comparisons.
#' 
#' The output is a dataframe containing estimated parameters and maximum
#' loglikelihood. The computed loglikelihood contains the factor q! m!/(q + m)!
#' where q is the number of species in the phylogeny and m is the number of
#' missing species, as explained in the supplementary material to Etienne et
#' al. 2012.
#' 
#' @param brtsM A set of branching times of the main clade in a phylogeny, all
#' positive
#' @param brtsS A set of branching times of the subclade in a phylogeny, all
#' positive
#' @param tsplit The branching time at which the lineage forming the subclade
#' branches off, positive
#' @param initparsopt The initial values of the parameters that must be
#' optimized
#' @param parsfix The values of the parameters that should not be optimized
#' @param idparsopt The ids of the parameters that must be optimized, e.g. 1:7
#' for all parameters.  The ids are defined as follows: \cr
#' id == 1 corresponds to lambda_M (speciation rate) of the main clade \cr
#' id == 2 corresponds to mu_M (extinction rate) of the main clade \cr
#' id == 3 corresponds to K_M (clade-level carrying capacity) of the main clade \cr id == 4 corresponds to
#' lambda_S (speciation rate) of the subclade \cr
#' id == 5 corresponds to mu_S (extinction rate) of the subclade \cr
#' id == 6 corresponds to t_d (the time of the key innovation)
#' @param idparsfix The ids of the parameters that should not be optimized,
#' e.g. c(1,3,4,6) if lambda and K should not be optimized, but only mu. In
#' that case idparsopt must be c(2,5,7). The default is to fix all parameters
#' not specified in idparsopt.
#' @param idparsnoshift The ids of the parameters that should not shift; This
#' can only apply to ids 4, 5 and 6, e.g. idparsnoshift = c(4,5) means that
#' lambda and mu have the same values before and after tshift
#' @param res sets the maximum number of species for which a probability must
#' be computed, must be larger than 1 + max(length(brtsM),length(brtsS))
#' @param ddmodel sets the model of diversity-dependence: \cr
#' \code{ddmodel == 1} : linear dependence in speciation rate with parameter K (= diversity
#' where speciation = extinction)\cr
#' \code{ddmodel == 1.3} : linear dependence in speciation rate with parameter K' (= diversity where speciation = 0)\cr
#' \code{ddmodel == 2} : exponential dependence in speciation rate with
#' parameter K (= diversity where speciation = extinction)\cr 
#' \code{ddmodel == 2.1} : variant of exponential dependence in speciation rate with offset at
#' infinity\cr
#' \code{ddmodel == 2.2} : 1/n dependence in speciation rate\cr
#' \code{ddmodel == 2.3} : exponential dependence in speciation rate with
#' parameter x (= exponent)\cr
#'\code{ddmodel == 3} : linear dependence in extinction rate \cr
#'\code{ddmodel == 4} : exponential dependence in
#' extinction rate \cr
#' \code{ddmodel == 4.1} : variant of exponential dependence in extinction rate
#' with offset at infinity \cr
#' \code{ddmodel == 4.2} : 1/n dependence in extinction rate with offset at
#' infinity \cr
#' @param missnumspec The number of species that are in the clade but missing
#' in the phylogeny. One can specify the sum of the missing species in main
#' clade and subclade or a vector c(missnumspec_M,missnumspec_S) with missing
#' species in main clade and subclade respectively.
#' @param cond Conditioning: \cr cond == 0 : no conditioning \cr cond == 1 :
#' conditioning on non-extinction of the phylogeny \cr
#' @param soc Sets whether stem or crown age should be used (1 or 2); stem age
#' only works when cond = 0
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
#' @param methode The method used in the ode solver, default is ode45
#' @param correction Sets whether the correction should be applied (TRUE) or
#' not (FALSE)
#' @param verbose Show the parameters and loglikelihood for every call to the
#' loglik function 
#' @return \item{lambda_M}{ gives the maximum likelihood estimate of lambda of
#' the main clade} \item{mu_M}{ gives the maximum likelihood estimate of mu of
#' the main clade} \item{K_M}{ gives the maximum likelihood estimate of K of
#' the main clade} \item{lambda_2}{ gives the maximum likelihood estimate of
#' lambda of the subclade} \item{mu_S}{ gives the maximum likelihood estimate
#' of mu of the subclade} \item{t_d}{ gives the time of the key innovation
#' event} \item{loglik}{ gives the maximum loglikelihood} \item{df}{ gives the
#' number of estimated parameters, i.e. degrees of feedom} \item{conv}{ gives a
#' message on convergence of optimization; conv = 0 means convergence}
#' @note The optimization may get trapped in local optima. Try different
#' starting values to search for the global optimum.
#' @author Rampal S. Etienne & Bart Haegeman
#' @seealso \code{\link{dd_MS_loglik}}, \code{\link{dd_ML}},
#' \code{\link{dd_KI_ML}}, \code{\link{dd_SR_ML}},
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#' 
#' cat("This will estimate parameters for two sets of branching times brtsM, brtsS\n")
#' cat("without conditioning.\n")
#' cat("The tolerance of the optimization is set high so runtime is fast in this example.\n")
#' cat("In real applications, use the default or more stringent settins for tol.\n")
#' brtsM = 4:10
#' brtsS = seq(0.1,3.5,0.7)
#' tsplit = 5
#' dd_MS_ML(brtsM = brtsM, brtsS = brtsS, tsplit = tsplit, idparsopt = c(1:3,6),
#'           initparsopt = c(0.885, 2e-14, 10, 4.001), idparsfix = NULL, parsfix = NULL,
#'           idparsnoshift = c(4,5), cond = 0, tol = c(3E-1,3E-1,3E-1))
#' 
#' @export dd_MS_ML
dd_MS_ML = function(brtsM,
    brtsS,
    tsplit,
    initparsopt = c(0.5, 0.1, 2 * (1 + length(brtsM) + length(brtsS) + sum(missnumspec)),(tsplit + max(brtsS))/2),
    parsfix = NULL,
    idparsopt = c(1:3,6),
    idparsfix = NULL,
    idparsnoshift = (1:6)[c(-idparsopt,(-1)^(length(idparsfix) != 0) * idparsfix)],
    res = 10*(1 + length(c(brtsM,brtsS)) + sum(missnumspec)),
    ddmodel = 1.3,
    missnumspec = 0,
    cond = 0,
    soc = 2,
    tol = c(1E-3, 1E-4, 1E-6),
    maxiter = 1000 * round((1.25)^length(idparsopt)),
    changeloglikifnoconv = FALSE,
    optimmethod = 'subplex',
    num_cycles = 1,
    methode = 'analytical',
    correction = FALSE,
    verbose = FALSE)
{
options(warn = -1)
brtsM = sort(abs(as.numeric(brtsM)),decreasing = TRUE)
brtsS = sort(abs(as.numeric(brtsS)),decreasing = TRUE)
if(cond == 1 & soc == 1)
{
   cat("Conditioning on survival of a clade with stem age currently not implemented.\n")
   out2 = data.frame(row.names = "results",lambda_M = -1, mu_M = -1, K = -1, lambda_S = -1, mu_S = -1, t_d = -1, loglik = -1, df = -1, conv = -1)
} else {
if(is.numeric(brtsM) == FALSE || is.numeric(brtsS) == FALSE)
{ 
   cat("The branching times should be numeric.\n")
   out2 = data.frame(row.names = "results",lambda_M = -1, mu_M = -1, K = -1, lambda_S = -1, mu_S = -1, t_d = -1, loglik = -1, df = -1, conv = -1)
} else {
idparsnoshift = sort(idparsnoshift)
idpars = sort(c(idparsopt,idparsfix,idparsnoshift))
if((prod(idpars == (1:6)) != 1) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)))
{
   cat("The parameters to be optimized, fixed and not shifted are incoherent.\n")
   out2 = data.frame(row.names = "results",lambda_M = -1, mu_M = -1, K = -1, lambda_S = -1, mu_S = -1, t_d = -1, loglik = -1, df = -1, conv = -1)
} else {
namepars = c("la_M","mu_M","K","la_S","mu_S","t_d")
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
pars2 = c(res,ddmodel,cond,tsplit,verbose,soc,correction,tol,maxiter)
optimpars = c(tol,maxiter)
initloglik = dd_MS_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsnoshift = idparsnoshift,pars2 = pars2,brtsM = brtsM,brtsS = brtsS,missnumspec = missnumspec, methode = methode)
cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
utils::flush.console()
if(initloglik == -Inf)
{
   cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
   out2 = data.frame(row.names = "results",lambda_M = -1, mu_M = -1, K = -1, lambda_S = -1, mu_S = -1, t_d = -1, loglik = -1, df = -1, conv = -1)
} else {
out = optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = dd_MS_loglik_choosepar,trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsnoshift = idparsnoshift,pars2 = pars2,brtsM = brtsM,brtsS = brtsS,missnumspec = missnumspec,methode = methode, num_cycles = num_cycles)
if(out$conv > 0)
{
   cat("Optimization has not converged. Try again with different initial values.\n")
   out2 = data.frame(row.names = "results",lambda_1 = -1, mu_1 = -1, K = -1, lambda_2 = -1, mu_2 = -1, t_d = -1, loglik = -1, df = -1, conv = unlist(out$conv))
} else {
MLtrpars = as.numeric(unlist(out$par))
MLpars = MLtrpars/(1-MLtrpars)
MLpars1 = rep(0,6)
MLpars1[idparsopt] = MLpars
if(length(idparsfix) != 0) {MLpars1[idparsfix] = parsfix }
if(length(idparsnoshift) != 0) { MLpars1[idparsnoshift] = MLpars1[idparsnoshift - 3] }
if(MLpars1[3] > 10^7){MLpars1[3] = Inf}
ML = as.numeric(unlist(out$fvalues))
out2 = data.frame(row.names = "results",lambda_M = MLpars1[1],mu_M = MLpars1[2],K = MLpars1[3], lambda_S = MLpars1[4], mu_S = MLpars1[5], t_d = MLpars1[6], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
if(out2$conv != 0 & changeloglikifnoconv == T) { out2$loglik = -Inf }
s1 = sprintf('Maximum likelihood parameter estimates: %f %f %f %f %f %f',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6])
s2 = sprintf('Maximum loglikelihood: %f',ML)
cat("\n",s1,"\n",s2,"\n")
out$par = list(MLpars1)
}
}
}
}
}
invisible(out2)
}

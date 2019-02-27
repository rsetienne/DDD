\name{dd_KI_loglik}
\alias{dd_KI_loglik}

\title{Loglikelihood for diversity-dependent diversification models with decoupling of a subclade from a main clade at time t = t_d }
\description{
  This function computes loglikelihood of a diversity-dependent diversification model
  for a given set of branching times and parameter values where the diversity-dependent dynamics
  of a subclade decouple from the dynamics of the main clade at time t_d, potentially
  accompanied by a shift in parameters.
}
\usage{dd_KI_loglik(pars1, pars2, brtsM, brtsS, missnumspec, methode = 'analytical')}
\arguments{
  \item{pars1}{Vector of parameters: \cr \cr
   \code{pars1[1]} corresponds to lambda_M (speciation rate) of the main clade \cr
   \code{pars1[2]} corresponds to mu_M (extinction rate) of the main clade \cr
   \code{pars1[3]} corresponds to K_M (clade-level carrying capacity) of the main clade \cr
   \code{pars1[4]} corresponds to lambda_M (speciation rate) of the subclade \cr
   \code{pars1[5]} corresponds to mu_S (extinction rate) of the subclade \cr
   \code{pars1[6]} corresponds to K_S (clade-level carrying capacity) of the subclade \cr
   \code{pars1[7]} corresponds to t_d (the time of decoupling)
}
  \item{pars2}{Vector of model settings: \cr \cr
    \code{pars2[1]} sets the maximum number of species for which a probability must be computed.
    This must be larger than 1 + missnumspec + length(brts). \cr \cr
    \code{pars2[2]} sets the model of diversity-dependence: \cr
    - \code{pars2[2] == 1} linear dependence in speciation rate with parameter K (= diversity where speciation = extinction)\cr
    - \code{pars2[2] == 1.3} linear dependence in speciation rate with parameter K' (= diversity where speciation = 0)\cr 
    - \code{pars2[2] == 2} exponential dependence in speciation rate with parameter K (= diversity where speciation = extinction)\cr
    - \code{pars2[2] == 2.1} variant of exponential dependence in speciation rate with offset at infinity\cr
    - \code{pars2[2] == 2.2} 1/n dependence in speciation rate\cr
    - \code{pars2[2] == 2.3} exponential dependence in speciation rate with parameter x (= exponent)\cr
    - \code{pars2[2] == 3} linear dependence in extinction rate \cr
    - \code{pars2[2] == 4} exponential dependence in extinction rate \cr
    - \code{pars2[2] == 4.1} variant of exponential dependence in extinction rate with offset at infinity\cr
    - \code{pars2[2] == 4.2} 1/n dependence in extinction rate\cr\cr
    \code{pars2[3]} sets the conditioning: \cr
    - \code{pars2[3] == 0} no conditioning \cr
    - \code{pars2[3] == 1} conditioning on non-extinction of the phylogeny \cr \cr
    \code{pars2[4]} sets the time of splitting of the branch that will decouple \cr \cr
    \code{pars2[5]} sets whether the parameters and likelihood should be shown on screen (1) or not (0) \cr \cr
    \code{pars2[6]} sets whether the first data point is stem age (1) or crown age (2) \cr \cr
    \code{pars2[7]} sets whether the old (incorrect) likelihood should be used (0), or whether new corrected version should be used (1)
}
  \item{brtsM}{A set of branching times of the main clade in the phylogeny, all positive }
  \item{brtsS}{A set of branching times of the subclade in the phylogeny, all positive }
  \item{missnumspec}{The number of species that are in the clade but missing in the phylogeny. One can specify the sum of the missing species in main clade and subclade or a vector c(missnumspec_M,missnumspec_S) with missing species in main clade and subclade respectively.}
  \item{methode}{The method used to solve the master equation, default is 'analytical' which uses matrix exponentiation; alternatively numerical ODE solvers can be used, such as 'lsoda' or 'ode45'. These were used in the package before version 3.1. } 
}
\value{The loglikelihood}
\references{
- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi: 10.1098/rspb.2011.1439 \cr
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi: 10.1086/667574
}
\author{Rampal S. Etienne & Bart Haegeman }
\seealso{
 \code{\link{dd_KI_ML}},
 \code{\link{dd_loglik}}
 \code{\link{dd_SR_loglik}}
}
\examples{
pars1 = c(0.25,0.12,25.51,1.0,0.16,8.61,9.8)
pars2 = c(200,1,0,18.8,1,2)
missnumspec = 0
brtsM = c(25.2,24.6,24.0,22.5,21.7,20.4,19.9,19.7,18.8,17.1,15.8,11.8,9.7,8.9,5.7,5.2)
brtsS = c(9.6,8.6,7.4,4.9,2.5)
dd_KI_loglik(pars1,pars2,brtsM,brtsS,missnumspec,methode = 'ode45')
}
\keyword{models}

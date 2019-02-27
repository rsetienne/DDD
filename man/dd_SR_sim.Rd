\name{dd_SR_sim}
\alias{dd_SR_sim}
\title{ Function to simulate the diversity-dependent diversification process with a shift in one or more of the parameters }
\description{ Simulating the diversity-dependent diversification process with a parameter shift at a certain time}
\usage{ dd_SR_sim(pars, age, ddmodel = 1) }
\arguments{
  \item{pars}{Vector of parameters: \cr \cr
    \code{pars[1]} corresponds to lambda1 (speciation rate before the rate shift) \cr
    \code{pars[2]} corresponds to mu1 (extinction rate before the rate shift)  \cr
    \code{pars[3]} corresponds to K1 (clade-level carrying capacity before the rate shift) \cr
    \code{pars[4]} corresponds to lambda2 (speciation rate after the rate shift) \cr
    \code{pars[5]} corresponds to mu2 (extinction rate after the rate shift)  \cr
    \code{pars[6]} corresponds to K2 (clade-level carrying capacity after the rate shift) \cr
    \code{pars[7]} corresponds to the time of shift 
  }
  \item{age}{ Sets the crown age for the simulation}
  \item{ddmodel}{ Sets the model of diversity-dependence: \cr
   \code{ddmodel == 1} : linear dependence in speciation rate with parameter K (= diversity where speciation = extinction)\cr
   \code{ddmodel == 1.3} : linear dependence in speciation rate with parameter K' (= diversity where speciation = 0)\cr
   \code{ddmodel == 2} : exponential dependence in speciation rate with parameter K (= diversity where speciation = extinction)\cr
   \code{ddmodel == 2.1} : variant of exponential dependence in speciation rate with offset at infinity\cr
   \code{ddmodel == 2.2} : 1/n dependence in speciation rate\cr
   \code{ddmodel == 2.3} : exponential dependence in speciation rate with parameter x (= exponent)\cr
   \code{ddmodel == 3} : linear dependence in extinction rate \cr
   \code{ddmodel == 4} : exponential dependence in extinction rate \cr
   \code{ddmodel == 4.1} : variant of exponential dependence in extinction rate with offset at infinity \cr
   \code{ddmodel == 4.2} : 1/n dependence in extinction rate with offset at infinity \cr
   \code{ddmodel == 5} : linear dependence in speciation and extinction rate
  }
}
\value{
  \item{ out }{ A list with the following four elements:
   The first element is the tree of extant species in phylo format \cr
   The second element is the tree of all species, including extinct species, in phylo format \cr
   The third element is a matrix of all species where \cr
   - the first column is the time at which a species is born \cr
   - the second column is the label of the parent of the species; positive and negative values only indicate whether the species belongs to the left or right crown lineage \cr
   - the third column is the label of the daughter species itself; positive and negative values only indicate whether the species belongs to the left or right crown lineage \cr
   - the fourth column is the time of extinction of the species If this equals -1, then the species is still extant.\cr
   The fourth element is the set of branching times of the tree of extant species.\cr
  }
}
\references{
- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi: 10.1098/rspb.2011.1439 \cr
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi: 10.1086/667574
}
\author{ Rampal S. Etienne }
\examples{ dd_SR_sim(c(0.2,0.1,20,0.2,0.1,40,5),10) }
\keyword{models}

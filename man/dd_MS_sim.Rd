\name{dd_MS_sim}
\alias{dd_MS_sim}
\title{ Function to simulate the macro-evolutionary succession process assuming diversity-dependent diversification }
\description{ Simulating a diversity-dependent diversification process where at a given time a new clade emerges with different inherent speciation rate and extinction rate }
\usage{ dd_MS_sim(pars, age, ddmodel = 1.3) }
\arguments{
  \item{pars}{Vector of parameters: \cr \cr
    \code{pars[1]} corresponds to lambda_M (speciation rate of the main clade) \cr
    \code{pars[2]} corresponds to mu_M (extinction rate of the main clade)  \cr
    \code{pars[3]} corresponds to K' (maximum number of species or a proxy for it in case of exponential decline in speciation rate) 
    \code{pars[4]} corresponds to lambda_S (speciation rate of the novel subclade) \cr
    \code{pars[5]} corresponds to mu_S (extinction rate)  \cr
    \code{pars[6]} tinn, the time the shift in rates occurs in the lineage leading to the subclade 
  }
  \item{age}{ Sets the crown age for the simulation }
  \item{ddmodel}{ Sets the model of diversity-dependence: \cr
   \code{ddmodel == 1.3} : linear dependence in speciation rate with parameter K' (= diversity where speciation = 0); ddmodel = 1 will be interpreted as this model \cr
   \code{ddmodel == 2.1} : variant of exponential dependence in speciation rate with offset at infinity; ddmodel = 2 will be interpreted as this model \cr
   \code{ddmodel == 2.2} : 1/n dependence in speciation rate\cr
   \code{ddmodel == 2.3} : exponential dependence in speciation rate with parameter x (= exponent)
  }
}
\value{
  \item{ out }{ A list with the following elements:
   The first element is the tree of extant species in phylo format \cr
   The second element is the tree of all species, including extinct species, in phylo format \cr
   The third element is a matrix of all species where \cr
   - the first column is the time at which a species is born \cr
   - the second column is the label of the parent of the species; positive and negative values only indicate whether the species belongs to the left or right crown lineage \cr
   - the third column is the label of the daughter species itself; positive and negative values only indicate whether the species belongs to the left or right crown lineage \cr
   - the fourth column is the time of extinction of the species \cr
   If the fourth element equals -1, then the species is still extant.\cr
   - the fifth column indicates whether the species belong to the main clade (0) or the subclade (1)\cr
   The fourth element is the subclade tree of extant species (without stem) \cr
   The fifth element is the subclade tree of all species (without stem) \cr
   The sixth element is the same as the first, except that it has attributed 0 for the main clade and 1 for the subclade\cr
   The seventh element is the same as the Second, except that it has attributed 0 for the main clade and 1 for the subclade\cr
   The sixth and seventh element will be NULL if the subclade does not exist (because it went extinct).
  }
}
\references{
- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi: 10.1098/rspb.2011.1439 \cr
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi: 10.1086/667574
}
\author{ Rampal S. Etienne }
\examples{ dd_MS_sim(c(0.2,0.1,20,0.1,0.05,4),10) }
\keyword{models}

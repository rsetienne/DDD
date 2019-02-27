\name{phylo2L}
\alias{phylo2L}
\title{ Function to convert phylogeny to a table with speciation and extinction events }
\description{ Converting a phylogeny to a table with speciation and extinction events }
\usage{ phylo2L(phy) }
\arguments{
  \item{ phy }{ A phylogeny of the phylo type }
}
\value{
  \item{L}{Matrix of events as produced by dd_sim: \cr \cr
   - the first column is the time at which a species is born in Mya\cr
   - the second column is the label of the parent of the species; positive and negative values indicate whether the species belongs to the left or right crown lineage \cr
   - the third column is the label of the daughter species itself; positive and negative values indicate whether the species belongs to the left or right crown lineage \cr
   - the fourth column is the time of extinction of the species; if the fourth element equals -1, then the species is still extant.
  }
}
\references{
- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi: 10.1098/rspb.2011.1439 \cr
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi: 10.1086/667574
}
\author{ Liang Xu }
\examples{
sim = dd_sim(c(0.2,0.1,20),10)
phy = sim$tas
L = phylo2L(phy)
phy2 = L2phylo(L, dropextinct = FALSE)
par(mfrow = c(1,3))
plot(phy)
plot(phy2)
plot(L2phylo(sim$L, dropextinct = FALSE))
}
\keyword{models}

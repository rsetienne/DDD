\name{roundn}
\alias{roundn}
\title{Rounds up in the usual manner }
\description{ The standard round function in R rounds x.5 to the nearest even integer. This is odd behavior that is corrected in roundn}
\usage{roundn(x, digits = 0)}
\arguments{
  \item{x}{Number to be rounded}
  \item{digits}{ Sets the number of decimals in rounding. }
}
\value{
  \item{n}{ A number }
}
\author{ Rampal S. Etienne }
\examples{
round(2.5)
roundn(2.5)
round(3.5)
roundn(3.5)
round(2.65,digits = 1)
roundn(2.65,digits = 1)
round(2.75,digits = 1)
roundn(2.75,digits = 1)
}
\keyword{models}

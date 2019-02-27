\name{sample2}
\alias{sample2}
\title{Takes samples in the usual manner }
\description{ The standard sample function in R samples from n numbers when x = n. This is unwanted behavior when the size of the vector to sample from changes dynamically. This is corrected in sample2}
\usage{sample2(x,size,replace = FALSE,prob = NULL) }
\arguments{
  \item{x}{A vector of one or more elements }
  \item{size}{A non-negative integer giving the number of items to choose. }
  \item{replace}{Should sampling be with replacement?}
  \item{prob}{A vector of probability weights for obtaining the elements of the vector being sampled.}
}
\value{
  \item{sam}{A vector of length \code{size} that is sampled from \code{x}. }
}
\author{Rampal S. Etienne }
\examples{
sample(x = 10,size = 5,replace = TRUE)
sample2(x = 10,size = 5,replace = TRUE)
}
\keyword{models}

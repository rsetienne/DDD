\name{optimizer}
\alias{optimizer}
\title{ Carries out optimization (finding a minimum) }
\description{ A wrapper to use several optimization routines, currently only 'simplex' (a method adopted from Matlab, or 'subplex', from the R package subplex). The function is called from several packages by the same author.}
\usage{
optimizer(optimmethod = 'simplex',
          optimpars = c(1E-4,1E-4,1E-6,1000),
          num_cycles = 1,
          fun,
          trparsopt,
          ...) }
\arguments{
  \item{optimmethod}{The method to use for optimization, either 'simplex' or 'subplex' }
  \item{optimpars}{Parameters of the optimization: relative tolerance in function arguments, relative tolerance in function value, absolute tolerance in function arguments, and maximum number of iterations }
  \item{num_cycles}{Number of cycles of the optimization. When set to Inf, the optimization will be repeated until the result is, within the tolerance, equal to the starting values, with a maximum of 5 cycles.}
  \item{fun}{Function to be optimized}
  \item{trparsopt}{Initial guess of the parameters to be optimized}
  \item{...}{Any other arguments of the function to be optimimzed, or settings of the optimization routine}
}
\value{
  \item{out}{ A list containing optimal function arguments (\code{par}, the optimal function value (\code{fvalues}) and whether the optimization converged (\code{conv})}.
}
\author{ Rampal S. Etienne }
\examples{
cat("No examples")
}
\keyword{models}

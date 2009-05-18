\name{NMDS}
\alias{NMDS}

\title{Nonmetric Multidimensional Scaling}
\description{
  Nonmetric multidimensional scaling by optimizing goodness of fit
  to a nonmetric hypothesis (Kruskal 1964)
}
\usage{
NMDS(D, k = 2, x, method = "BFGS", ...)
}

\arguments{
  \item{D}{ Dissimilarity structure. }
  \item{k}{ Number of dimensions. }
  \item{x}{ Optional starting configuration. A random start is used
   if missing.}
  \item{method}{Minimization method for \code{\link{optim}}. }
  \item{\dots}{ Other arguments passed to \code{\link{optim}}. }
}
\details{
  Function consists of stress function, derivative function of the stress,
  and call to \code{\link{optim}}. 
}
\value{
  Function returns a result object of \code{\link{optim}}. For most items, 
  see the documentation of \code{\link{optim}}. The following items have
  a specific meaning in \code{NMDS}:
  \item{par}{The final configuration.}
  \item{value}{The final stress.}
}
\references{ 
  Kruskal, J. B. (1964) Nonmetric multidimensional scaling: a numerical
  method. \emph{Psychometrika} 29, 115--129.
 }
\author{ Jari Oksanen }

\seealso{ \code{\link[MASS]{isoMDS}} function (\pkg{MASS} package) 
  implements the same algorithm but written in C. See \code{\link{optim}}
  for the optimization routine. }

\keyword{ multivariate }


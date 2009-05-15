\name{CAwa1}
\Rdversion{1.1}
\alias{CAwa1}

\title{
  Correspondence Analysis.
}
\description{
  First axis of Correspondence Analysis using weighted averaging, also
  known as Reciprocal Averaging (Hill 1973).
}
\usage{
CAwa1(Y)
}

\arguments{
  \item{Y}{
    Data matrix.
}
}

\details{ This is a poor algorithm to perform CA, and is included only
  to show how simple it is. Correspondence Analysis was originally
  introduced in ecology as a result of this algorithm, and called
  Reciprocal Averaging (Hill 1973).  The function only finds one axis.
  It would be possible to find later axes by orthogonalizing against
  previous axes, but if several axes are needed, it is better to use
  other algorithms from the beginning.  }
\value{
  \item{eig}{Eigenvalue of the first axis.}
  \item{u}{Row scores of the first axis.}
  \item{v}{Column scores of the first axis.}
}

\references{
  Hill, M.O. (1973) Reciprocal averaging: an eigenvector method of
  ordination. \emph{J. Ecol.} 61, 237--249.
}
\author{
Jari Oksanen.
}

\seealso{
 \code{\link{CA}} and \code{\link{CAeig}} are better algorihtms for
 CA. \code{\link{PCApot1}} is a similar poor algorithm for PCA.
}

\keyword{ multivariate }


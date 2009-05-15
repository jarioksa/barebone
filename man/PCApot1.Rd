\name{PCApot1}
\Rdversion{1.1}
\alias{PCApot1}

\title{
  Principal Component Analysis

}
\description{
  First axis of Principal Component Analysis using power method.
  }
\usage{
PCApot1(Y, scale = FALSE)
}

\arguments{
  \item{Y}{
    Data matrix.
}
  \item{scale}{
    Scale columns to unit variance.
}
}
\details{
  This is the poorest algorithm of performing PCA, and is included only
  to show how simple it is. The function only finds one axis. It would
  be possible to find later axes by orthogonalizing against previous
  axes, but that is not worthwhile with an algorithm as poor as this
  one. 
}
\value{
  The function returns an object with items:
  \item{eig}{Eigenvalue of the first axis.}
  \item{ueig}{Rowscores scaled by eigenvalues.}
  \item{v}{Orthonormal column scores.}
}
\references{
?
}
\author{
  Jari Oksanen.
}
\seealso{
  \code{\link{PCA}} and \code{\link{PCAeig}} for better algorithms.
}

\keyword{ multivariate }

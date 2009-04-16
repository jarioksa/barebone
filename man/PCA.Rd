\name{PCA}
\alias{PCA}

\title{ Principal Components Analysis }
\description{
  Principal Components analysis via singular value decomposition.
}
\usage{
PCA(Y, scale = FALSE)
}

\arguments{
  \item{Y}{ Data Matrix }
  \item{scale}{Scale columns to unit variance.}
}



\value{
 Function returns an object from \code{\link{svd}} with items:
  \item{d}{Singular values which are the square roots of eigenvalues.}
  \item{u}{Orthonormal row scores.}
  \item{v}{Orthonormal column scores.}
}
\author{ Jari Oksanen }

\note{ 
  Functions \code{\link{prcomp}} and \code{\link[vegan]{rda}}
  (package \pkg{vegan}) use a similar algorithm.  
}

\seealso{ \code{\link{PCAeig}}, \code{\link{svd}}, \code{\link{prcomp}}. }

\keyword{ multivariate }


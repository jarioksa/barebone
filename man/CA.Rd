\name{CA}
\alias{CA}

\title{Correspondence Analysis }
\description{
 Correspondence Analysis via singular value decomposition.
}
\usage{
CA(Y)
}

\arguments{
  \item{Y}{ Data Matrix }
}

\value{
 Function returns an object from \code{\link{svd}} with items:
  \item{d}{Singular values which are the square roots of eigenvalues.}
  \item{u}{Orthonormal row scores.}
  \item{v}{Orthonormal column scores.}
}

\references{ 
Greenacre, M. J. (1984). Theory and applications of correspondence
     analysis. Academic Press, London.
}
\author{ Jari Oksanen }

\note{
Function \code{\link[vegan]{cca}} (\pkg{vegan}) package uses similar algorithm.
}
\seealso{ \code{\link{pCCA}}, \code{\link{svd}}, 
\code{\link[vegan]{cca}}. }

\keyword{ multivariate }


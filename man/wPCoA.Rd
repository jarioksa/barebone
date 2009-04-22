\name{wPCoA}
\alias{wPCoA}

\title{Weighted Principal Coordinate Analysis }
\description{
 Function performs Principal Coordinate Analysis with row weights.
}
\usage{
wPCoA(D, w)
}
\arguments{
  \item{D}{ Dissimilarity or distance structure as in 
    \code{\link{PCoA}}. }
  \item{w}{ Weights for observations. }
}

\value{
 Function returns an object from \code{\link{eigen}}.
  \item{values}{Eigenvalues}
  \item{vectors}{Orthonormal row scores.}
}

\author{ Jari Oksanen }

\seealso{ \code{\link{PCoA}}, \code{\link[vegan]{wcmdscale}}. }

\keyword{ multivariate }


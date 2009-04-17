\name{CAeig}
\alias{CAeig}

\title{ Correspondence Analysis }
\description{
 Correspondence Analysis via eigen decomposition of symmetric
 crossproduct matrix.
}
\usage{
CAeig(Y)
}

\arguments{
  \item{Y}{ Data matrix. }
}

\value{
  The function returns an object from \code{\link{eigen}} enhanced
  with row scores.
  \item{values}{Eigenvalues.}
  \item{vectors}{Weighted orthonormal column scores.}
  \item{u}{Weighted rthonormal row scores.}

}
\references{
Nishisato, S. 1980. Analysis of categorical data: Dual scaling and its
application. Univ Toronto Press.  
}
\author{Jari Oksanen}

\seealso{\code{\link{CA}}, \code{\link{PCAeig}}, \code{\link{eigen}}. }

\keyword{ multivariate }


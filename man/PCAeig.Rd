\name{PCAeig}
\alias{PCAeig}

\title{ Principal Components Analysis }
\description{
 Principal Components Analysis via eigen decomposition of crossprodcut matrix.
}
\usage{
PCAeig(Y, scale = FALSE)
}

\arguments{
  \item{Y}{ Data matrix. }
  \item{scale}{Scale columns to unit variance. }
}

\value{
  The function returns an object from \code{\link{eigen}} enhanced
  with row scores.
  \item{values}{Eigenvalues.}
  \item{vectors}{Orthonormal colum scores.}
  \item{u}{Orthonormal row scores.}
}

\author{ Jari Oksanen }
\note{ 
 Function \code{\link{princomp}} uses a similar algorithm.
}
\seealso{ \code{\link{PCA}}, \code{\link{eigen}}, \code{\link{princomp}}. }

\keyword{ multivariate}


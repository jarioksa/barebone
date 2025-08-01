% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCAeig.R
\name{PCAeig}
\alias{PCAeig}
\title{Principal Components Analysis}
\usage{
PCAeig(Y, scale = FALSE)
}
\arguments{
\item{Y}{Data matrix.}

\item{scale}{Scale columns to unit variance.}
}
\value{
The function returns an object from \code{\link{eigen}} enhanced
with row scores.  \item{values}{Eigenvalues.} \item{vectors}{Orthonormal
colum scores.} \item{u}{Orthonormal row scores.}
}
\description{
Principal Components Analysis via eigen decomposition of crossprodcut
matrix.
}
\note{
Function \code{\link{princomp}} uses a similar algorithm.
}
\seealso{
\code{\link{PCA}}, \code{\link{eigen}}, \code{\link{princomp}}.
}
\author{
Jari Oksanen
}
\keyword{multivariate}

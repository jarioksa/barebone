\name{pCCA}
\alias{pCCA}

\title{ [Partial] [Constrained] Correspondence Analysis. }
\description{
 [Partial] [Constrained] Correspondence Analysis via singular 
  value and QR decomposition.
}
\usage{
pCCA(Y, X = NULL, Z = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ Dependent data Matrix. }
  \item{X}{ Model matrix of constraints (can be missing). }
  \item{Z}{ Model matrix of conditions (can be missing). }
}

\value{
  later\dots
}

\author{ Jari Oksanen}
\note{
  Function \code{\link[vegan]{cca}} (\pkg{vegan} package) is similar. 
}
\seealso{ \code{\link{pRDA}}, \code{\link{svd}}, \code{\link{qr}},
 \code{\link[vegan]{cca}}. }

\keyword{ multivariate }


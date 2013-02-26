\name{wRDA}
\alias{wRDA}

\title{  Weighted [Partial] Redundancy Analysis and Principal Components Analysis }
\description{
  Weighted [Partial] Redundancy Analysis via QR and singular value
  decomposition or Principal Components Analysis via singular value
  decomposition.
}
\usage{
wRDA(Y, X = NULL, Z = NULL, scale = FALSE, w = NULL)
}

\arguments{
  \item{Y}{ Dependent Data Matrix. }
  \item{X}{ Model matrix of constraints (can be missing). }
  \item{Z}{ Model matrix of constraints (can be missing). }
  \item{scale}{ Scale dependent matrix. }
  \item{w}{Row weights.}
}

\value{
later\dots
}

\author{ Jari Oksanen }
\note{ 
  Function \code{\link{pRDA}} is a similmar function for non-weighted
  [partial] RDA.
}

\seealso{ 
  \code{\link{pRDA}}, \code{\link{pCCA}}, \code{\link[vegan]{rda}}.
}

\keyword{ multivariate }


\name{pRDA}
\alias{pRDA}
\alias{pRDAalt}

\title{ [Partial] Redundancy Analysis and Principal Components Analysis }
\description{
  [Partial] Redundancy Analysis via QR and singular value decomposition
  or Principal Components Analysis via singular value decomposition.
}
\usage{
pRDA(Y, X = NULL, Z = NULL, scale = FALSE)
pRDAalt(Y, X = NULL, Z = NULL, scale = FALSE)
}

\arguments{
  \item{Y}{ Dependent Data Matrix. }
  \item{X}{ Model matrix of constraints (can be missing). }
  \item{Z}{ Model matrix of constraints (can be missing). }
  \item{scale}{ Scale dependent matrix. }
}

\value{
later\dots
}

\author{ Jari Oksanen }
\note{ Function \code{\link[vegan]{rda}} is essentially similar.
  Function \code{\link{RDA}} is a simpler basic function that only
  performs Redundancy Analysis. Function \code{pRDAalt} takes explicit
  residuals of \code{X} regressed on \code{Z} (similarly as CANOCO for
  Windows software), whereas \code{pRDA} trusts QR decomposition to
  handle the orthogonality.}

\seealso{ \code{\link{RDA}}, \code{\link{pCCA}}, \code{\link{svd}},
  \code{\link{qr}}, \code{\link[vegan]{rda}}. }

\keyword{ multivariate }


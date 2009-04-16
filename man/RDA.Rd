\name{RDA}
\alias{RDA}

\title{ Redundancy Analysis }
\description{
  Redundancy Analysis via QR and singular value decompositions.
}
\usage{
RDA(Y, X, scale = FALSE)
}

\arguments{
  \item{Y}{ Dependent data matrix. }
  \item{X}{ Model matrix of constraints. }
  \item{scale}{ Scaling of dependent data. }
}

\value{
  later\dots
}

\author{ Jari Oksanen}
\note{ 
  Function \code{\link{pRDA}} is a more versatile an complicated
  variant, that also can do partial RDA or PCA as special cases.
}
\seealso{ \code{\link{pRDA}}, \code{\link{svd}}, \code{\link{qr}},
  \code{\link[vegan]{rda}}. }

\keyword{ multivariate }


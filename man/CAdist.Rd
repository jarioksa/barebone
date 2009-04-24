\name{CAdist}
\alias{CAdist}

\title{ Correspondence Analysis }
\description{
  Correspondence analysis as a weighted Principal Coordinate Analysis
  of Chi-squared distances.
}
\usage{
CAdist(Y)
}

\arguments{
  \item{Y}{ Data matrix. }
}

\value{
 The function returns an object from \code{\link{wPCoA}} which
 in turn returns an \code{\link{eigen}} result amended with
 species scores.
  \item{values }{Eigenvalues.}
  \item{vectors }{Weighted orthonormal row scores.}
  \item{v}{Weighted orthonormal column scores.}
}

\author{ Jari Oksanen }

\keyword{ multivariate }


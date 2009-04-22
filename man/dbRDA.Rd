\name{dbRDA}
\alias{dbRDA}

\title{Distance-based Redundancy Analysis}
\description{
  Distance-based Redundancy Analysis is a Redundancy Analysis on
  the eigenvalue-scaled result of Principal Coordinates Analysis.
}
\usage{
dbRDA(D, X, Z)
}

\arguments{
  \item{D}{ Dissimilarity structure of dependent data 
    as in \code{\link{PCoA}}. }
  \item{X}{ Model matrix of constraints (can be missing).  }
  \item{Z}{ Model matrix of conditions (can be missing). }
}

\value{
  The function returns the same object as \code{\link{dbRDA}}.
}

\author{ Jari Oksanen }

\seealso{ \code{\link{PCoA}}, \code{\link{pRDA}}, \code{\link[vegan]{capscale}}.  }

\keyword{multivariate}


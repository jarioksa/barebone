% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbRDA.R
\name{dbRDA}
\alias{dbRDA}
\title{Distance-based Redundancy Analysis}
\usage{
dbRDA(D, X)
}
\arguments{
\item{D}{Dissimilarity structure of dependent data as in \code{\link{PCoA}}.}

\item{X}{Model matrix of constraints.}
}
\value{
The function returns the same object as \code{\link{RDA}}.
}
\description{
Distance-based Redundancy Analysis is a Redundancy Analysis of
Gower-standardized dissimilarity matrix.
}
\references{
McArdle, B.H. & Anderson, M.J. (2001).  Fitting
    multivariate models to community data: a comment on
    distance-based redundancy analysis.  _Ecology_ 82, 290-297.
}
\author{
Jari Oksanen

#' @seealso \code{\link{PCoA}}, \code{\link{RDA}},
\code{\link[vegan]{dbrda}}.
}
\keyword{multivariate}

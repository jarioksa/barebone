\name{uht2beta}
\alias{uht2beta}
\title{ Gaussian response coefficients }

\description{ The function converts parameters of a Gaussian response curve
  into polynomial linear coefficients.  }

\usage{
uht2beta(u, h, t)
}
\arguments{
  \item{u}{ Optimum (location) of the curve. }
  \item{h}{ Expected height of the curve at optimum. }
  \item{t}{ Tolerance (standard deviation of the normal density). }
}

\details{ This function is useful for generating Gaussian
  species responses to environmental gradients using e.g.
  the \code{\link{rspecies}} function.
}

\value{
  A 3 element vector for regression coefficients:
  first is the intecept, second is for the linear term,
  third is for the squared term.
}
\author{ Peter Solymos }
\examples{
uht2beta(25, 10, 20)
uht2beta(75, 20, 5)
}
\keyword{ datagen }
\keyword{ multivariate }

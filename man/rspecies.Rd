\name{rspecies}
\alias{rspecies}
\title{ Community Matrix Generation }
\description{
  The function generates a random community matrix, where species' abundances and environmental relationships can be specified. Abundances are modelled by individual linear models with log link and the Poisson distribution.
}
\usage{
rspecies(n, spp, b = rep(1, spp), x = rep(1, n))
}
\arguments{
  \item{n}{ The number of sites (rows) in the data matrix. }
  \item{spp}{ The number of species (columns) in the data matrix. }
  \item{b}{ Coefficients to be used ofr determining species' environmental relationships. It can be a vector of length \code{spp} if only intercepts are provided, or a list with \code{spp} elements (one for each species, each element must be of the same length \eqn{p} for the same covariates, including the intercept). }
  \item{x}{ A vector (if only intercepts are provided) or a model matrix with \code{n} rows and \eqn{p} columns.  }
}
\details{
  The linear model is simply \code{log.lambda = x %*% b}, and then the returned vector for one species is \code{rpois(n, exp(log.lambda))}, following a standard Poisson GLM.

  The model matrix should be in the usual format, providing design variables for factors. The formula interface of the \link{\code{model.matrix}} function can be applied (see Examples).
}
\value{
  A community data matrix with \code{n} rows and \code{spp} columns, containing counts.
}
\author{ Peter Solymos }
\examples{
## 5 sites and 2 species
rspecies(5,2, c(2,4))
## intercept and 1 continuous covariate
x <- rnorm(10)
rspecies(10,2, list(c(2,1), c(1,2)), model.matrix(~x))
## intercept and 1 fixed factor
x <- rep(0:1, each = 5)
rspecies(10,2, list(c(1,-1), c(1,1)), model.matrix(~x))
if (require(vegan) {
## number of sites and species
n <- 1000
spp <- 50
## parameters for lognormal RAD
mu <- 5    # mean for rlnorm
sigma <- 1 # sd for rlnorm
x1 <- rspecies(n, spp, rnorm(spp, mu, sigma))
## parameters for the exponential and broken stick RADs
alpha <- 0.1 # decay
J <- 100     # community size
x2 <- rspecies(n, spp, log(J * alpha * (1 - alpha)^(1:spp - 1)))
x3 <- rspecies(n, spp, log(J * cumsum(1/spp:1)/spp))
## fitting RADs
(rf1 <- radfit(colSums(x1))) # lognormal
(rf2 <- radfit(colSums(x2))) # exponential
(rf3 <- radfit(colSums(x3))) # brokenstick
## plot the results
opar <- par(mfrow=c(1,3))
plot(rf1)
plot(rf2)
plot(rf3)
par(opar)
}
}
\keyword{ datagen }
\keyword{ multivariate }

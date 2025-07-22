#' Nonmetric Multidimensional Scaling
#'
#' Nonmetric multidimensional scaling by optimizing goodness of fit to a
#' nonmetric hypothesis (Kruskal 1964)
#'
#' Function consists of stress function, derivative function of the stress, and
#' call to \code{\link{optim}}.
#'
#' @param D Dissimilarity structure.
#' @param k Number of dimensions.
#' @param u Optional starting configuration. A random start is used if missing.
#' @param method Minimization method for \code{\link{optim}}.
#' @param \dots Other arguments passed to \code{\link{optim}}.
#' @return Function returns a result object of \code{\link{optim}}. For most
#' items, see the documentation of \code{\link{optim}}. The following items
#' have a specific meaning in \code{NMDS}: \item{par}{The final configuration.}
#' \item{value}{The final stress.}
#' @author Jari Oksanen
#' @seealso \code{\link[MASS]{isoMDS}} function (\pkg{MASS} package) implements
#' the same algorithm but written in C. See \code{\link{optim}} for the
#' optimization routine and \code{\link{isoreg}} for isotonic regression.
#' @references Kruskal, J. B. (1964) Nonmetric multidimensional scaling: a
#' numerical method. \emph{Psychometrika} 29, 115--129.
#' @keywords multivariate
#'
#' @importFrom stats isoreg dist runif optim
#'
#' @export NMDS
NMDS <-
function(D, k = 2, u, method = "L-BFGS-B", ...)
{
    stress <- function(u, D, k) {
        u <- matrix(u, ncol=k)
        y <- dist(u)
        ord <- order(D, y)
        s <- isoreg(y[ord])
        sqrt(sum((s$y - s$yf)^2)/sum(s$y^2))
    }
    dstress <- function(u, D, k) {
        u <- matrix(u, ncol=k)
        dx <- u
        y <- dist(u)
        ord <- order(D, y)
        s <- isoreg(y[ord])
        yf <- s$yf[order(ord)]
        Sstar <- sum((y-yf)^2)
        Tstar <- sum(y^2)
        S <- sqrt(Sstar/Tstar)
        dmat <- (S/Sstar*(y-yf) - S/Tstar*y)/y
        dmat <-  as.matrix(dmat)
        for (l in 1:k)
            dx[,l] <- rowSums(dmat * outer(u[,l], u[,l], "-"))
        as.vector(dx)
    }
    D <- round(D, 12)
    if (missing(u))
        u <- runif(k * attr(D, "Size"))
    out <- optim(u, stress, dstress, D = D, k = k, method = method, ...)
    out$par <- matrix(out$par, ncol = k)
    out
}


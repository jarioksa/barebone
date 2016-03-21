#' Principal Component Analysis
#' 
#' First axis of Principal Component Analysis using power method.
#' 
#' This is the poorest algorithm of performing PCA, and is included only to
#' show how simple it is. The function only finds one axis. It would be
#' possible to find later axes by orthogonalizing against previous axes, but
#' that is not worthwhile with an algorithm as poor as this one.
#' 
#' @param Y Data matrix.
#' @param scale Scale columns to unit variance.
#' @return The function returns an object with items: \item{eig}{Eigenvalue of
#' the first axis.} \item{ueig}{Rowscores scaled by eigenvalues.}
#' \item{v}{Orthonormal column scores.}
#' @author Jari Oksanen.
#' @seealso \code{\link{PCA}} and \code{\link{PCAeig}} for better algorithms.
#' @references ?
#' @keywords multivariate
#'
#' @importFrom stats runif
#' 
#' @export PCApot1
`PCApot1` <- function(Y, scale = FALSE)
{
    EPS <- sqrt(.Machine$double.eps)
    Y <- as.matrix(scale(Y, center = TRUE, scale = scale))
    v <- runif(ncol(Y))
    v <- v - mean(v)
    v <- v/sqrt(sum(v^2))
    eig <- 0
    repeat {
        u <- Y %*% v
        v <- t(Y) %*% u
        ss <- sqrt(sum(v^2))
        v <- v/ss
        if (abs(ss-eig) < EPS)
            break
        eig <- ss
    }
    list(eig = eig, ueigen  = u, v = v)
}

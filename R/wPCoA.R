#' Weighted Principal Coordinate Analysis
#' 
#' Function performs Principal Coordinate Analysis with row weights.
#' 
#' 
#' @param D Dissimilarity or distance structure as in \code{\link{PCoA}}.
#' @param w Weights for observations.
#' @return Function returns an object from \code{\link{eigen}}.
#' \item{values}{Eigenvalues} \item{vectors}{Orthonormal row scores.}
#' @author Jari Oksanen
#' @seealso \code{\link{PCoA}}, \code{\link[vegan]{wcmdscale}}.
#' @keywords multivariate
#'
#' @importFrom stats weighted.mean
#' 
#' @export wPCoA
`wPCoA` <-
    function(D, w)
{
    M <- as.matrix(D^2)
    cnt <- apply(M, 2, weighted.mean, w = w)
    M <- scale(M, center = cnt, scale = FALSE)
    cnt <- apply(M, 1, weighted.mean, w = w)
    M <- t(scale(t(M), center = cnt, scale = FALSE))
    M <- diag(sqrt(w)) %*% M %*% diag(sqrt(w))
    SOL <- eigen(-M/2, symmetric = TRUE)
    SOL$vectors <- diag(1/sqrt(w)) %*% SOL$vectors
    SOL
}

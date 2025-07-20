#' Principal Coordinates Analysis
#'
#' Principal Coordinates Analysis (PCoA) is also known as metric or classic
#' multidimensional scaling.
#'
#' Function does not scale the result vectors by eigenvalues, and returns
#' all non-zero eigenvalues, including negative eigenvalues, and associated
#' axes.
#'
#' @param D Dissimilarity or distance structure. This either should inherit
#' from \code{\link{dist}} or have a similar \code{\link{as.matrix}} method
#' producing a symmetric matrix.
#' @return Function returns an object from \code{\link{eigen}}.
#' \item{values}{Eigenvalues} \item{vectors}{Orthonormal row scores.}
#' @author Jari Oksanen
#' @seealso \code{\link{cmdscale}}, \code{\link{PCAeig}}.
#' @references Gower, J. C. (1966) Some distance properties of latent root and
#' vector methods used in multivariate analysis.  \emph{Biometrika} \bold{53},
#' 325--328.
#'
#' Mardia, K. V., Kent, J. T. and Bibby, J. M. (1979).  Chapter 14 of
#' \emph{Multivariate Analysis}, London: Academic Press.
#' @keywords multivariate
#' @export PCoA
`PCoA` <-
    function(D)
{
    M <- as.matrix(D^2)
    M <- scale(M, center = TRUE, scale = FALSE)
    M <- t(scale(t(M), center = TRUE, scale = FALSE))
    SOL <- eigen(-M/2, symmetric = TRUE)
    nonzero <- abs(SOL$values) > sqrt(.Machine$double.eps)
    SOL$values <- SOL$values[nonzero]
    SOL$vectors <- SOL$vectors[, nonzero]
    SOL
}

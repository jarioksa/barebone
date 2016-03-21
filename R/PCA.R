#' Principal Components Analysis
#' 
#' Principal Components analysis via singular value decomposition.
#' 
#' 
#' @param Y Data Matrix
#' @param scale Scale columns to unit variance.
#' @return Function returns an object from \code{\link{svd}} with items:
#' \item{d}{Singular values which are the square roots of eigenvalues.}
#' \item{u}{Orthonormal row scores.} \item{v}{Orthonormal column scores.}
#' @note Functions \code{\link{prcomp}} and \code{\link[vegan]{rda}} (package
#' \pkg{vegan}) use a similar algorithm.
#' @author Jari Oksanen
#' @seealso \code{\link{PCAeig}}, \code{\link{svd}}, \code{\link{prcomp}}.
#' @keywords multivariate
#' @export PCA
`PCA` <-
    function(Y, scale = FALSE)
{
    Ybar <- scale(Y, center = TRUE, scale = scale)
    svd(Ybar)
}

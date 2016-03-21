#' Principal Components Analysis
#' 
#' Principal Components Analysis via eigen decomposition of crossprodcut
#' matrix.
#' 
#' 
#' @param Y Data matrix.
#' @param scale Scale columns to unit variance.
#' @return The function returns an object from \code{\link{eigen}} enhanced
#' with row scores.  \item{values}{Eigenvalues.} \item{vectors}{Orthonormal
#' colum scores.} \item{u}{Orthonormal row scores.}
#' @note Function \code{\link{princomp}} uses a similar algorithm.
#' @author Jari Oksanen
#' @seealso \code{\link{PCA}}, \code{\link{eigen}}, \code{\link{princomp}}.
#' @keywords multivariate
#' @export PCAeig
`PCAeig` <-
    function(Y, scale = FALSE)
{
    Y <- scale(as.matrix(Y), center = TRUE, scale = scale)
    YY <- crossprod(Y)
    SOL <- eigen(YY, symmetric = TRUE)
    SOL$u <- Y %*% SOL$vectors %*% diag(1/sqrt(SOL$values))
    SOL
}

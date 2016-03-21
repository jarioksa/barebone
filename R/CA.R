#' Correspondence Analysis
#' 
#' Correspondence Analysis via singular value decomposition.
#' 
#' 
#' @param Y Data Matrix
#' @return Function returns an object from \code{\link{svd}} with items:
#' \item{d}{Singular values which are the square roots of eigenvalues.}
#' \item{u}{Orthonormal row scores.} \item{v}{Orthonormal column scores.}
#' @note Function \code{\link[vegan]{cca}} (\pkg{vegan}) package uses similar
#' algorithm.
#' @author Jari Oksanen
#' @seealso \code{\link{pCCA}}, \code{\link{svd}}, \code{\link[vegan]{cca}}.
#' @references Greenacre, M. J. (1984). Theory and applications of
#' correspondence analysis. Academic Press, London.
#' @keywords multivariate
#'
#' @export CA
`CA` <-
    function(Y)
{
    Y <- Y/sum(Y)
    r <- rowSums(Y)
    c <- colSums(Y)
    Ybar <- (Y - r %o% c)/sqrt(r %o% c)
    SOL <- svd(Ybar)
    SOL$u <- diag(1/sqrt(r)) %*% SOL$u
    SOL$v <- diag(1/sqrt(c)) %*% SOL$v
    SOL
}

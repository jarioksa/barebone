#' Correspondence Analysis.
#' 
#' First axis of Correspondence Analysis using weighted averaging, also known
#' as Reciprocal Averaging (Hill 1973).
#' 
#' This is a poor algorithm to perform CA, and is included only to show how
#' simple it is. Correspondence Analysis was originally introduced in ecology
#' as a result of this algorithm, and called Reciprocal Averaging (Hill 1973).
#' The function only finds one axis.  It would be possible to find later axes
#' by orthogonalizing against previous axes, but if several axes are needed, it
#' is better to use other algorithms from the beginning.
#' 
#' @param Y Data matrix.
#' @return \item{eig}{Eigenvalue of the first axis.} \item{u}{Row scores of the
#' first axis.} \item{v}{Column scores of the first axis.}
#' @author Jari Oksanen.
#' @seealso \code{\link{CA}} and \code{\link{CAeig}} are better algorihtms for
#' CA. \code{\link{PCApot1}} is a similar poor algorithm for PCA.
#' @references Hill, M.O. (1973) Reciprocal averaging: an eigenvector method of
#' ordination. \emph{J. Ecol.} 61, 237--249.
#' @keywords multivariate
#'
#' @importFrom vegan wascores
#' @importFrom stats rnorm
#' 
#' @export CAwa1
`CAwa1` <-
    function(Y)
{
    EPS <- sqrt(.Machine$double.eps)
    u <- rnorm(nrow(Y))
    eig <- 0
    repeat {
        v <- wascores(u, Y, expand = TRUE)
        v <- sweep(v, 1, attr(v, "centre"))
        u <- wascores(v, t(Y), expand = TRUE)
        if (abs(attr(u, "shrinkage") - eig) < EPS)
            break
        eig <- attr(u, "shrinkage")
    }
    list(eig = attr(u, "shrinkage"), u = u, v = v)
}

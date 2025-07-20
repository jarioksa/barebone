#' Constrained Correspondence Analysis.
#' 
#' First axis of Constrained Correspondence Analysis using weighted averaging.
#' 
#' This is a poor algorithm to perform CCA, and is included only to show how
#' simple it is. Constrained Correspondence Analysis was originally introduced
#' in ecology as a result of this algorithm (ter Braak 1986). Palmer (1993)
#' gives a very lucid explanation of the algorithm and its two kind of row
#' scores, known as Weighted Averages (WA) scores and Linear Combination (LC)
#' scores.  The function only finds one axis.  It would be possible to find
#' later axes by orthogonalizing against previous axes, but if several axes are
#' needed, it is better to use other algorithms from the beginning.
#' 
#' @param Y Data matrix.
#' @param X Model matrix of constraints.
#' @return \item{eig}{Eigenvalue of the first axis.} \item{w}{WA scores for
#' rows of the first axis.} \item{u}{LC scores for rows of the first axis.}
#' \item{v}{Column scores of the first axis.}
#' @author Jari Oksanen.
#' @seealso \code{\link{CA}} and \code{\link{CAeig}} are better algorihtms for
#' CA. \code{\link{PCApot1}} is a similar poor algorithm for PCA.
#' @references Palmer, M. W. (1993) Putting things in even better order: The
#' advantages of canonical correspondence analysis. \emph{Ecology} 74,
#' 2215--2230.
#' 
#' ter Braak, C. J. F. (1986) Canonical correspondence analysis: a new
#' eigenvector technique of multivariate direct gradient analysis.
#' \emph{Ecology} 67, 1167--1179.
#' @keywords multivariate
#'
#' @importFrom stats rnorm fitted lm
#' @importFrom vegan wascores
#' 
#' @export CCAwa1
`CCAwa1` <-
    function(Y, X)
{
    EPS <- sqrt(.Machine$double.eps)
    X <- as.matrix(X)
    wts <- rowSums(Y)
    u <- rnorm(nrow(Y))
    eig <- 0
    repeat {
        v <- wascores(u, Y, expand = TRUE)
        w <- wascores(v, t(Y), expand = TRUE)
        u <- fitted(lm(w ~  X, weights = wts))
        if (abs(attr(v, "shrinkage") - eig) < EPS)
            break
        eig <- attr(v, "shrinkage")
    }
    list(eig = attr(v, "shrinkage"), w = w, u = u, v = v)
}

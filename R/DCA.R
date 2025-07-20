#' Detrended Correspondence Analysis
#'
#' Correspondence Analysis where axes are detrended by
#' \code{\link{loess}}.
#'
#' Function performes detrended correspondence analysis but it differs
#' essentially from the classic \code{\link[vegan]{decorana}}. Not
#' only is the detrending method different, but \code{decorana} makes
#' much more than detrending. Most importantly it attempts to rescale
#' axes to constant response width (see \pkg{vegan} function
#' \code{\link[vegan]{tolerance}}) which often has a stronger effect
#' than detrending.
#'
#' If you want to inspect the steps of \code{decorana}, function
#' \code{rdecorana} in github package \pkg{natto} implements the
#' function in plain \R{} which is much easier to follow than the
#' original Fortran code.
#'
#' @param Y Data matrix.
#' @param pairwise Detrend axis pairwise against each previous axis,
#'     or against all previous axes together.
#' @param \dots Other arguments passed to \code{\link{loess}}.
#' @author Jari Oksanen
#' @seealso \code{\link[vegan]{decorana}}.
#' @return Function returns a list with items:
#' \item{DCA.eig}{Eigenvalue as estimated during detrending
#'     steps.}
#' \item{eig}{Eigenvalue as estimated by shrinking of
#'     scores in final solution, similarly as in \code{\link{CA}}.}
#' \item{u, v}{Row and column scores.}
#' @keywords multivariate
#' @importFrom vegan wascores eigengrad
#' @importFrom graphics points
#' @importFrom stats runif weighted.mean residuals loess fitted
#' @export
`DCA` <-
    function(Y, pairwise = FALSE, ...)
{
    EPS <- sqrt(.Machine$double.eps)
    CYCLES <- 200
    Y <- Y/sum(Y)
    r <- rowSums(Y)
    c <- colSums(Y)
    V <- matrix(0, ncol(Y), 4)
    U <- matrix(0, nrow(Y), 4)
    EIG <- numeric(4)
    v <- runif(ncol(Y))
    eig <- eig.old <- 0
    for (k in 1:4) {
        cycles <- 0
        repeat {
            v <- v - weighted.mean(v, w=c)
            v <- v/sqrt(sum(c*v*v))
            u <- wascores(v, t(Y))
            if (k > 1) {
                if (pairwise && k > 2) {
                   for(kk in c(1:(k-1), (k-2):1))
                       u <- residuals(loess(u ~ U[,kk], weights = r, ...))
                }
                else
                    u <- residuals(lo <- loess(u ~ U[,1:(k-1)], weights = r,
                                           normalize = FALSE,  ...))
                plot(u ~ U[,1], cex=0.3)
                points(U[,1], fitted(lo), pch=16, col=2)
            }
            vprime <- wascores(u, Y)
            eig <- abs(sum(c*v*vprime))
            if ((tol <- abs(eig - eig.old)) < EPS || CYCLES <
                (cycles <- cycles+1)) {
                break
            } else {
                v <- vprime
                eig.old <- eig
            }
        }
        EIG[k] <- eig
        V[,k] <- v
        U[,k] <- wascores(v, t(Y))
        v <- runif(length(v))
        if (tol > EPS)
            warning("residual ", formatC(tol, digits=4), " larger than tolerance in axis ", k)
    }
    eig0 <- eigengrad(V, t(Y))
    list(DCA.eig = EIG, eig = eig0, u = U, v = V)
}

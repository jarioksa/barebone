#' Weighted [Partial] Redundancy Analysis and Principal Components Analysis
#' 
#' Weighted [Partial] Redundancy Analysis via QR and singular value
#' decomposition or Principal Components Analysis via singular value
#' decomposition.
#' 
#' 
#' @param Y Dependent Data Matrix.
#' @param X Model matrix of constraints (can be missing).
#' @param Z Model matrix of constraints (can be missing).
#' @param scale Scale dependent matrix.
#' @param w Row weights.
#' @return later\dots{}
#' @note Function \code{\link{pRDA}} is a similmar function for non-weighted
#' [partial] RDA.
#' @author Jari Oksanen
#' @seealso \code{\link{pRDA}}, \code{\link{pCCA}}, \code{\link[vegan]{rda}}.
#' @keywords multivariate
#' @export wRDA
`wRDA` <-
    function(Y, X = NULL, Z = NULL, scale = FALSE, w = NULL, cw = NULL)
{
    if (is.null(w))
        w <- rep(1/nrow(Y), nrow(Y))
    else
        w <- w/sum(w)
    if (is.null(cw))
        cw <- rep(1/ncol(Y), ncol(Y))
    else
        cw <- cw/sum(cw)
    tot <- sum(Y)/(nrow(Y) - 1)
    wscale <- function(x, scale = FALSE, w = NULL) {
        wc <- apply(x, 2, weighted.mean, w = w)
        x <- sweep(x, 2, wc, "-")
        if (scale) {
            wsd <- sqrt(colSums(w * x^2) / (nrow(Y) - 1))
            x <- sweep(x, 2, wsd, "/")
        }
        sweep(x, 1, sqrt(w), "*")
    }
    Y <- wscale(as.matrix(Y), scale = scale, w = w)
    Y <- sweep(Y, 2, sqrt(cw), "*")
    Y <- sqrt(tot) * Y
    if (!is.null(Z)) {
        Z <- wscale(as.matrix(Z), scale = FALSE, w = w)
        QZ <- qr(Z)
        Y <- qr.resid(QZ, Y)
    }
    if (!is.null(X)) {
        X <- wscale(as.matrix(X), scale = FALSE, w = w)
        X <- cbind(X, Z)
        Q <- qr(X)
        RDA <- svd(qr.fitted(Q, Y))
        RDA$d <- RDA$d^2
        RDA$u <- sqrt(1/w) * RDA$u
        RDA$v <- sqrt(1/cw) * RDA$v
        RDA$w <- Y %*% RDA$v %*% diag(1/RDA$d)
        Y <- qr.resid(Q, Y)
    } else {
        RDA <- NULL
    }
    RES <- svd(Y)
    RES$d <- RES$d^2
    RES$u <- sqrt(1/w) * RES$u
    RES$v <- sqrt(1/cw) * RES$v
    list(RDA = RDA, RES = RES)
}


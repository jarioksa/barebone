`wRDA` <-
    function(Y, X = NULL, Z = NULL, scale = FALSE, w = NULL)
{
    if (is.null(w))
        w <- rep(1/nrow(Y), nrow(Y))
    else
        w <- w/sum(w)
    wscale <- function(x, scale = FALSE, w = NULL) {
        wc <- apply(x, 2, weighted.mean, w = w)
        x <- sweep(x, 2, wc, "-")
        if (scale) {
            wsd <- sqrt(colSums(w * x^2) / (nrow(Y) - 1))
            x <- sweep(x, 2, wsd, "/")
        }
        sweep(x, 1, w, "*")
    }
    Y <- wscale(as.matrix(Y), scale = scale, w = w)
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
        RDA$w <- Y %*% RDA$v %*% diag(1/RDA$d)
        Y <- qr.resid(Q, Y)
    } else {
        RDA <- NULL
    }
    RES <- svd(Y)
    list(RDA = RDA, RES = RES)
}


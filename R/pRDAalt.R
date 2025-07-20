#' @rdname pRDA
#' @export
`pRDAalt` <-
function(Y, X = NULL, Z = NULL, scale = FALSE)
{
    Y <- scale(as.matrix(Y), center = TRUE, scale = scale)
    if (!is.null(Z)) {
        Z <- scale(as.matrix(Z), center = TRUE, scale = FALSE)
        QZ <- qr(Z)
        Y <- qr.resid(QZ, Y)
    }
    if (!is.null(X)) {
        X <- scale(as.matrix(X), center = TRUE, scale = FALSE)
        if (!is.null(Z))
            X <- qr.resid(QZ, X)
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


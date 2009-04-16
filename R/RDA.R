`RDA` <-
function(Y, X, scale = FALSE)
{
    Y <- scale(as.matrix(Y), center = TRUE, scale = scale)
    X <- scale(as.matrix(X), center = TRUE, scale = FALSE)
    Q <- qr(X)
    RDA <- svd(qr.fitted(Q, Y))
    RDA$w <- Y %*% RDA$v %*% diag(1/RDA$d)
    RES <- svd(qr.resid(Q, Y))
    list(RDA = RDA, Resid = RES)
}


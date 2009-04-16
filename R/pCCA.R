`pCCA` <-
    function(Y, X = NULL, Z = NULL)
{
    Y <- as.matrix(Y)/sum(Y)
    r <- rowSums(Y)
    c <- colSums(Y)
    Y <- (Y - r %o% c)/sqrt(r %o% c)
    if (!is.null(Z)) {
        Z <- as.matrix(Z)
        rmean <- apply(Z, 2, weighted.mean, w = r)
        Z <- scale(Z, center = rmean, scale = FALSE)
        Z <- diag(sqrt(r)) %*% Z
        Y <- qr.resid(qr(Z), Y)
    }
    if (!is.null(X)) {
        X <- as.matrix(X)
        rmean <- apply(X, 2, weighted.mean, w = r)
        X <- scale(X, center = rmean, scale = FALSE)
        X <- diag(sqrt(r)) %*% X
        X <- cbind(X, Z)
        Q <- qr(X)
        CCA <- svd(qr.fitted(Q, Y))
        CCA$w <- Y %*% CCA$v
        CCA$w <- diag(1/sqrt(r)) %*% CCA$w %*% diag(1/CCA$d)
        CCA$u <- diag(1/sqrt(r)) %*% CCA$u
        CCA$v <- diag(1/sqrt(c)) %*% CCA$v
        Y <- qr.resid(Q, Y)
    } else {
        CCA <- NULL
    }
    RES <- svd(Y)
    RES$u <- diag(1/sqrt(r)) %*% RES$u
    RES$v <- diag(1/sqrt(c)) %*% RES$v
    list(CCA = CCA, RES = RES)
}
    

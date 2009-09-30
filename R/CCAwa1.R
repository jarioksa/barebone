`CCAwa1` <-
    function(Y, X)
{
    require(vegan) || stop("requires wascores function in vegan")
    EPS <- sqrt(.Machine$double.eps)
    X <- as.matrix(X)
    u <- rnorm(nrow(Y))
    eig <- 0
    repeat {
        v <- wascores(u, Y, expand = TRUE)
        w <- wascores(v, t(Y), expand = TRUE)
        u <- fitted(lm(w ~  X))
        if (abs(attr(v, "shrinkage") - eig) < EPS)
            break
        eig <- attr(v, "shrinkage")
    }
    list(eig = attr(v, "shrinkage"), w = w, u = u, v = v)
}

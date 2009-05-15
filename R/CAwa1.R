`CAwa1` <-
    function(Y)
{
    require(vegan) || stop("requires wascores function in vegan")
    EPS <- sqrt(.Machine$double.eps)
    u <- rnorm(nrow(Y))
    eig <- 0
    repeat {
        v <- wascores(u, Y, expand = TRUE)
        u <- wascores(v, t(Y), expand = TRUE)
        if (abs(attr(u, "shrinkage") - eig) < EPS)
            break
        eig <- attr(u, "shrinkage")
    }
    list(eig = attr(u, "shrinkage"), u = u, v = v)
}

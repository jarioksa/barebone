`PCApot1` <- function(Y, scale = FALSE)
{
    EPS <- sqrt(.Machine$double.eps)
    Y <- as.matrix(scale(Y, center = TRUE, scale = scale))
    v <- runif(ncol(Y))
    v <- v - mean(v)
    v <- v/sqrt(sum(v^2))
    eig <- 0
    repeat {
        u <- Y %*% v
        v <- t(Y) %*% u
        ss <- sqrt(sum(v^2))
        v <- v/ss
        if (abs(ss-eig) < EPS)
            break
        eig <- ss
    }
    list(eig = eig, ueigen  = u, v = v)
}

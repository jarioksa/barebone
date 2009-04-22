`wPCoA` <-
    function(D, w)
{
    M <- as.matrix(D^2)
    cnt <- apply(M, 2, weighted.mean, w = w)
    M <- scale(M, center = cnt, scale = FALSE)
    cnt <- apply(M, 1, weighted.mean, w = w)
    M <- t(scale(t(M), center = cnt, scale = FALSE))
    M <- diag(sqrt(w)) %*% M %*% diag(sqrt(w))
    SOL <- eigen(-M/2, symmetric = TRUE)
    SOL$vectors <- diag(1/sqrt(w)) %*% SOL$vectors
    SOL
}

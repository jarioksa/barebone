`CAdist` <-
    function(Y)
{
    Y <- as.matrix(Y)/sum(Y)
    r <- rowSums(Y)
    c <- colSums(Y)
    Ybar <- diag(1/r) %*% Y %*% diag(1/sqrt(c)) 
    D <- dist(Ybar)
    SOL <- wPCoA(D, r)
    SOL$v <- diag(1/c) %*% t(Y) %*% SOL$vectors %*% diag(1/sqrt(SOL$values))
    ## The next gives the same result with Ybar instead of Y
    ##SOL$v <- diag(1/sqrt(c)) %*% t(Ybar) %*% diag(r) %*% SOL$vectors %*% diag(1/sqrt(SOL$values))
    SOL
}

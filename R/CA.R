`CA` <-
    function(Y)
{
    Y <- Y/sum(Y)
    r <- rowSums(Y)
    c <- colSums(Y)
    Y <- (Y - r %o% c)/sqrt(r %o% c)
    SOL <- svd(Y)
    SOL$u <- diag(1/sqrt(r)) %*% SOL$u
    SOL$v <- diag(1/sqrt(c)) %*% SOL$v
    SOL
}

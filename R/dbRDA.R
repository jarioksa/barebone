`dbRDA` <-
    function(D, X, Z)
{
    SOL <- PCoA(D)
    pos <- SOL$values > sqrt(.Machine$double.eps)
    Lambda <- SOL$values[pos]
    Y <- SOL$vectors[, pos] %*% diag(sqrt(Lambda))
    pRDA(Y, X, Z)
}

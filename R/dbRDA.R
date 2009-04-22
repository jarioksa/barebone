`dbRDA` <-
    function(D, X, Z)
{
    SOL <- PCoA(D)
    npos <- SOL$values > sqrt(.Machine$double.eps)
    Lambda <- SOL$values[npos]
    Y <- SOL$vectors[, npos] %*% diag(sqrt(Lambda))
    pRDA(Y, X, Z)
}

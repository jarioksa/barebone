`PCoA` <-
    function(D)
{
    M <- as.matrix(D^2)
    M <- scale(M, center = TRUE, scale = FALSE)
    M <- t(scale(t(M), center = TRUE, scale = FALSE))
    SOL <- eigen(-M/2, symmetric = TRUE)
    nonzero <- abs(SOL$values) > sqrt(.Machine$double.eps)
    SOL$values <- SOL$values[nonzero]
    SOL$vectors <- SOL$vectors[, nonzero]
    SOL
}

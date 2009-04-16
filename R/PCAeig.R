`PCAeig` <-
    function(Y, scale = FALSE)
{
    Y <- scale(as.matrix(Y), center = TRUE, scale = scale)
    YY <- crossprod(Y)
    SOL <- eigen(YY, symmetric = TRUE)
    SOL$u <- Y %*% SOL$vectors %*% diag(1/sqrt(SOL$values))
    SOL
}

`CAeig` <-
function (Y) 
{
    Y <- as.matrix(Y)/sum(Y)
    r <- rowSums(Y)
    c <- colSums(Y)
    YY <- diag(1/sqrt(c)) %*% t(Y) %*% diag(1/r) %*% Y %*% diag(1/sqrt(c))
    SOL <- eigen(YY, symmetric = TRUE)
    SOL$vectors <- diag(1/sqrt(c)) %*% SOL$vectors
    SOL$u <- Y %*% SOL$vectors
    SOL$u <- diag(1/r) %*% SOL$u %*% diag(1/sqrt(SOL$values))
    SOL$u <- SOL$u[,-1]
    SOL$vectors <- SOL$vectors[,-1]
    SOL$values <- SOL$values[-1]
    SOL
}


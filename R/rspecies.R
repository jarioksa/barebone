rspecies <-
function(n, spp, b=rep(1, spp), x=rep(1, n))
{
    log.lambda <- lapply(1:spp, function(z) data.matrix(x) %*% b[[z]])
    y <- sapply(log.lambda, function(z) rpois(n, exp(z)))
    y
}

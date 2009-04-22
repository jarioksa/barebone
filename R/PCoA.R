`PCoA` <-
    function(D)
{
    M <- as.matrix(D^2)
    M <- scale(M, center = TRUE, scale = FALSE)
    M <- t(scale(t(M), center = TRUE, scale = FALSE))
    SOL <- eigen(-M/2, symmetric = TRUE)
    SOL
}

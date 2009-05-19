NMDS <-
function(D, k = 2, x, method = "BFGS", ...)
{
    stress <- function(x, D, k) {
        x <- matrix(x, ncol=k)
        y <- dist(x)
        ord <- order(D, y)
        s <- isoreg(y[ord])
        sqrt(sum((s$y - s$yf)^2)/sum(s$y^2))
    }
    dstress <- function(x, D, k) {
        x <- matrix(x, ncol=k)
        dx <- x
        y <- dist(x)
        ord <- order(D, y)
        s <- isoreg(y[ord])
        yf <- s$yf[order(ord)]
        Sstar <- sum((y-yf)^2)
        Tstar <- sum(y^2)
        S <- sqrt(Sstar/Tstar)
        dmat <- (S/Sstar*(y-yf) - S/Tstar*y)/y
        dmat <-  as.matrix(dmat) 
        for (l in 1:k)
            dx[,l] <- rowSums(dmat * outer(x[,l], x[,l], "-"))
        as.vector(dx)
    }
    if (missing(x))
        x <- runif(k * attr(D, "Size"))
    out <- optim(x, stress, dstress, D = D, k = k, method = method, ...)
    out$par <- matrix(out$par, ncol = k)
    out
}


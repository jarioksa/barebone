NMDS <-
function(D, k = 2, u, method = "BFGS", ...)
{
    stress <- function(u, D, k) {
        u <- matrix(u, ncol=k)
        y <- dist(u)
        ord <- order(D, y)
        s <- isoreg(y[ord])
        sqrt(sum((s$y - s$yf)^2)/sum(s$y^2))
    }
    dstress <- function(u, D, k) {
        u <- matrix(u, ncol=k)
        dx <- u
        y <- dist(u)
        ord <- order(D, y)
        s <- isoreg(y[ord])
        yf <- s$yf[order(ord)]
        Sstar <- sum((y-yf)^2)
        Tstar <- sum(y^2)
        S <- sqrt(Sstar/Tstar)
        dmat <- (S/Sstar*(y-yf) - S/Tstar*y)/y
        dmat <-  as.matrix(dmat) 
        for (l in 1:k)
            dx[,l] <- rowSums(dmat * outer(u[,l], u[,l], "-"))
        as.vector(dx)
    }
    D <- round(D, 12)
    if (missing(u))
        u <- runif(k * attr(D, "Size"))
    out <- optim(u, stress, dstress, D = D, k = k, method = method, ...)
    out$par <- matrix(out$par, ncol = k)
    out
}


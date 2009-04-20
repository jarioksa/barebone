rspecies <-
function(n, spp, b=rep(1, spp), x=rep(1, n))
{
    log.lambda <- lapply(1:spp, function(z) data.matrix(x) %*% b[[z]])
    y <- sapply(log.lambda, function(z) rpois(n, exp(z)))
    y
}

## brokenstick
## cumsum(1/spp:1)/spp

## exponential
## alpha * (1 - alpha)^(1:spp - 1)

## lognormal
## mu * sigma^rnorm(spp)



n <- 1000
spp <- 50
mu <- 5
sigma <- 1
alpha <- 0.1
J <- 100

#x1 <- rspecies(n, spp, mu * sigma^rnorm(spp))
x1 <- rspecies(n, spp, rnorm(spp, 5, 1))
x2 <- rspecies(n, spp, exp(J * alpha * (1 - alpha)^(1:spp - 1)))
x3 <- rspecies(n, spp, exp(J * cumsum(1/spp:1)/spp))

(rf1 <- radfit(colSums(x1)))
(rf2 <- radfit(colSums(x2)))
(rf3 <- radfit(colSums(x3)))

opar <- par(mfrow=c(1,3))
plot(rf1)
plot(rf2)
plot(rf3)
par(opar)

explore <-
function(x, FUN, i, j, k=1)
{
    rout <- list(list(list()))
    for (ii in 1:length(i)) {
        for (jj in 1:length(j)) {
            for (kk in 1:k)
                rout[[ii]][[jj]][[kk]] <- FUN(x, i[ii], j[jj])
        }
    }
    rval
}

fill <- function(x) sum(x > 0) / prod(dim(x))
uht2beta <-
function(u, h, t)
{
    b2 <- -1 / (2 * t^2)
    b1 <- -2 * u * b2
    b0 <- log(h) + b1^2 / (4 * b2)
    c(b0, b1, b2)
}

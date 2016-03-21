#' Redundancy Analysis
#' 
#' Redundancy Analysis via QR and singular value decompositions.
#' 
#' 
#' @param Y Dependent data matrix.
#' @param X Model matrix of constraints.
#' @param scale Scaling of dependent data.
#' @return later\dots{}
#' @note Function \code{\link{pRDA}} is a more versatile an complicated
#' variant, that also can do partial RDA or PCA as special cases.
#' @author Jari Oksanen
#' @seealso \code{\link{pRDA}}, \code{\link{svd}}, \code{\link{qr}},
#' \code{\link[vegan]{rda}}.
#' @keywords multivariate
#' @export RDA
`RDA` <-
function(Y, X, scale = FALSE)
{
    Y <- scale(as.matrix(Y), center = TRUE, scale = scale)
    X <- scale(as.matrix(X), center = TRUE, scale = FALSE)
    Q <- qr(X)
    RDA <- svd(qr.fitted(Q, Y))
    RDA$w <- Y %*% RDA$v %*% diag(1/RDA$d)
    RES <- svd(qr.resid(Q, Y))
    list(RDA = RDA, Resid = RES)
}


#' Distance-based Redundancy Analysis
#'
#' Distance-based Redundancy Analysis is a Redundancy Analysis of
#' Gower-standardized dissimilarity matrix.
#'
#' @param D Dissimilarity structure of dependent data as in \code{\link{PCoA}}.
#' @param X Model matrix of constraints.
#' @return The function returns the same object as \code{\link{RDA}}.
#' @author Jari Oksanen
#'
#' #' @seealso \code{\link{PCoA}}, \code{\link{RDA}},
#' \code{\link[vegan]{dbrda}}.
#' @references McArdle, B.H. & Anderson, M.J. (2001).  Fitting
#'     multivariate models to community data: a comment on
#'     distance-based redundancy analysis.  _Ecology_ 82, 290-297.
#' @keywords multivariate
#' @export dbRDA
`dbRDA` <-
    function(D, X)
{
    M <- as.matrix(D^2)
    M <- scale(M, scale = FALSE)
    M <- t(scale(t(M), scale = FALSE))
    M <- -M/2
    Q <- qr(scale(X, scale=FALSE))
    fit <- qr.fitted(Q, M)
    fit <- qr.fitted(Q, t(fit))
    RDA <- eigen(fit, symmetric = TRUE)
    pos <- RDA$values > sqrt(.Machine$double.eps)
    RDA$values <- RDA$values[pos]
    RDA$vectors <- RDA$vectors[, pos]
    RDA$w <- M %*% RDA$vectors %*% diag(1/RDA$values)
    M <- qr.resid(Q, M)
    M <- qr.resid(Q, t(M))
    RES <- eigen(M, symmetric = TRUE)
    nz <- abs(RES$value) > sqrt(.Machine$double.eps)
    RES$values <- RES$values[nz]
    RES$vectors <- RES$vectors[, nz]
    list(RDA = RDA, Resid = RES)
}

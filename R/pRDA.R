#' [Partial] Redundancy Analysis and Principal Components Analysis
#'
#' [Partial] Redundancy Analysis via QR decompositon and singular
#' value decomposition or Principal Components Analysis via singular
#' value decomposition.
#'
#'
#' @aliases pRDA pRDAalt
#' @param Y Dependent Data Matrix.
#' @param X Model matrix of constraints (can be missing).
#' @param Z Model matrix of constraints (can be missing).
#' @param scale Scale dependent matrix.
#' @return Function returns the a list of components \code{RDA} and
#'     \code{RES} from constrained and residual unconstrained
#'     analysis. Thease are the output of \code{\link{svd}}, but
#'     component \code{RDA} is amended with WA scores \code{w}.
#' @note Function \code{\link[vegan]{rda}} (\pkg{vegan}) is
#'     essentially similar.  Function \code{\link{RDA}} is a simpler
#'     basic function that only performs Redundancy Analysis. Function
#'     \code{pRDAalt} takes explicit residuals of \code{X} regressed
#'     on \code{Z} (similarly as CANOCO for Windows software), whereas
#'     \code{pRDA} trusts QR decomposition to handle the
#'     orthogonality.
#' @author Jari Oksanen
#' @seealso \code{\link{RDA}}, \code{\link{pCCA}}, \code{\link{svd}},
#'     \code{\link{qr}}, \code{\link[vegan]{rda}}.
#' @keywords multivariate
#' @export pRDA
`pRDA` <-
function(Y, X = NULL, Z = NULL, scale = FALSE)
{
    Y <- scale(as.matrix(Y), center = TRUE, scale = scale)
    if (!is.null(Z)) {
        Z <- scale(as.matrix(Z), center = TRUE, scale = FALSE)
        QZ <- qr(Z)
        Y <- qr.resid(QZ, Y)
    }
    if (!is.null(X)) {
        X <- scale(as.matrix(X), center = TRUE, scale = FALSE)
        X <- cbind(X, Z)
        Q <- qr(X)
        RDA <- svd(qr.fitted(Q, Y))
        RDA$w <- Y %*% RDA$v %*% diag(1/RDA$d)
        Y <- qr.resid(Q, Y)
    } else {
        RDA <- NULL
    }
    RES <- svd(Y)
    list(RDA = RDA, RES = RES)
}


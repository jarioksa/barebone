#' Distance-based Redundancy Analysis
#'
#' Distance-based Redundancy Analysis is a Redundancy Analysis on the
#' eigenvalue-scaled result of Principal Coordinates Analysis.
#'
#'
#' @param D Dissimilarity structure of dependent data as in \code{\link{PCoA}}.
#' @param X Model matrix of constraints (can be missing).
#' @param Z Model matrix of conditions (can be missing).
#' @return The function returns the same object as \code{\link{pRDA}}.
#' @author Jari Oksanen
#'
#' @note This function does not actually perform distance-based RDA,
#'     but a simple variant called \code{\link[vegan]{capscale}} in
#'     \pkg{vegan}. The function will be fixed ASAP.
#'
#' @seealso \code{\link{PCoA}}, \code{\link{pRDA}},
#' \code{\link[vegan]{capscale}}.
#' @references Legendre, P. & Anderson, M.J. (1999). Canonical analysis of
#' principal coordinates: a useful method of constrained ordination for
#' ecology. \emph{Ecology} 84, 511--525.
#' @keywords multivariate
#' @export dbRDA
`dbRDA` <-
    function(D, X = NULL, Z = NULL)
{
    SOL <- PCoA(D)
    pos <- SOL$values > sqrt(.Machine$double.eps)
    Lambda <- SOL$values[pos]
    Y <- SOL$vectors[, pos] %*% diag(sqrt(Lambda))
    pRDA(Y, X, Z)
}

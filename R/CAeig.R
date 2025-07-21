#' Correspondence Analysis
#'
#' Correspondence Analysis via eigen decomposition of symmetric crossproduct
#' matrix.
#'
#'
#' @param Y Data matrix.
#' @return The function returns an object from \code{\link{eigen}} enhanced
#' with row scores.  \item{values}{Eigenvalues.} \item{vectors}{Weighted
#' orthonormal column scores.} \item{u}{Weighted orthonormal row scores.}
#' @author Jari Oksanen
#' @seealso \code{\link{CA}}, \code{\link{PCAeig}}, \code{\link{eigen}}.
#' @references Nishisato, S. 1980. Analysis of categorical data: Dual scaling
#' and its application. Univ Toronto Press.
#' @keywords multivariate
#' @export CAeig
`CAeig` <-
function (Y)
{
    Y <- as.matrix(Y)/sum(Y)
    r <- rowSums(Y)
    c <- colSums(Y)
    YY <- diag(1/sqrt(c)) %*% t(Y) %*% diag(1/r) %*% Y %*% diag(1/sqrt(c))
    SOL <- eigen(YY, symmetric = TRUE)
    SOL$vectors <- diag(1/sqrt(c)) %*% SOL$vectors
    SOL$u <- Y %*% SOL$vectors
    SOL$u <- diag(1/r) %*% SOL$u %*% diag(1/sqrt(SOL$values))
    SOL$u <- SOL$u[,-1]
    SOL$vectors <- SOL$vectors[,-1]
    SOL$values <- SOL$values[-1]
    SOL$values[SOL$values < 0] <- 0
    SOL
}


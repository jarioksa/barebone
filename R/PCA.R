`PCA` <-
    function(Y, scale = FALSE)
{
    Y <- scale(Y, center = TRUE, scale = scale)
    svd(Y)
}

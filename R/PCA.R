`PCA` <-
    function(Y, scale = FALSE)
{
    Ybar <- scale(Y, center = TRUE, scale = scale)
    svd(Ybar)
}

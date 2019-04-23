distEuclidean <- function (x, centers) 
{
  if (ncol(x) != ncol(centers)) 
    stop(sQuote("x"), " and ", sQuote("centers"), " must have the same number of columns")
  z <- matrix(0, nrow = nrow(x), ncol = nrow(centers))
  for (k in 1:nrow(centers)) {
    z[, k] <- sqrt(colSums((t(x) - centers[k, ])^2))
  }
  z
}
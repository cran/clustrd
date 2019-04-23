fitted.tuneclus <- function(object, mth = c("centers", "classes"), ...) 
{
  mth <- match.arg(mth)
  object$clusobjbest$centroid = data.frame(object$clusobjbest$centroid)
  if (mth == "centers") 
    object$clusobjbest$centroid[object$clusobjbest$cluster, , drop = FALSE]
  else object$clusobjbest$cluster
}
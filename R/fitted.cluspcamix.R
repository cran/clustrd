fitted.cluspcamix <- function(object, mth = c("centers", "classes"), ...) 
{
  mth <- match.arg(mth)
  
  if (mth == "centers") 
    if (object$centroid != 0) {
      object$centroid[object$cluster, , drop = FALSE]
    } else {
      object$centroid[object$cluster]
    }
  else { 
    if (object$centroid != 0) {
      object$cluster
    } else {
      object$centroid[object$cluster]
    }
  }
}
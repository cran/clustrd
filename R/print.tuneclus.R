print.tuneclus <- function(x, ...) {
  size = x$clusobjbest$size
  csize = round((table(x$clusobjbest$cluster)/sum(table(x$clusobjbest$cluster)))*100,digits=1)
  k = x$nclusbest
  d = x$ndimbest
  x$clusobjbest$centroid = data.frame(x$clusobjbest$centroid)
  cluspca = FALSE
  try(if (x$clusobjbest$center) { cluspca = TRUE }, silent = TRUE)
  
  if (cluspca == TRUE) {  
    if (x$clusobjbest$center == TRUE) {
      centering = "mean centered" } else {
        centering = "not centered"
      }
    if (x$clusobjbest$scale == TRUE) {
      scaling = "standardized" } else {
        scaling = "unstandardized"
      }
    
    tt = paste('(',csize,'%)',sep="")
    cs = paste(size, tt, sep = " ", collapse = ", ")
    cat(paste("\nThe best solution was obtained for ",k ," clusters of sizes ", paste(cs, collapse = ", ")," in ",d ," dimensions, for a cluster quality criterion value of ",round(x$critbest,3), ". Variables were ", centering, " and ", scaling,".", "\n", sep = ""))
  } else {
    
    tt = paste('(',csize,'%)',sep="")
    cs = paste(size, tt, sep = " ", collapse = ", ")
    cat(paste("\nThe best solution was obtained for ",k ," clusters of sizes ", paste(cs, collapse = ", ")," in ",d ," dimensions, for a cluster quality criterion value of ",round(x$critbest,3), ".", "\n", sep = ""))
  }
  
  cat("\nCluster quality criterion values across the specified range of clusters (rows) and dimensions (columns):\n")
  print(x$critgrid)
  
  cat("\nCluster centroids:\n")
  xcent = data.frame(round(x$clusobjbest$centroid,4))
  for (i in 1:k) {
    rownames(xcent)[i] = paste("Cluster",i)
  }
  for (i in 1:ncol(xcent)) {
    colnames(xcent)[i] = paste0("Dim.",i)
  }
  print(xcent)
 # attc = data.frame(round(x$clusobjbest$attcoord,4))
#  cat("\nVariable scores:\n")
 # for (i in 1:ncol(attc)) {
#    colnames(attc)[i] = paste0("Dim.",i)
#  }
#  print(attc)
  
  cat("\nWithin cluster sum of squares by cluster:\n")
  #resid <- x$obscoord - fitted(x) 
  #tot.withinss <- ss(resid)
  #print(tot.withinss)
  betweenss <- ss(x$clusobjbest$centroid[x$clusobjbest$cluster,]) # or 
  #betweenss <- ss(fitted(x))
  withinss <- sapply(split(as.data.frame(x$clusobjbest$obscoord), x$clusobjbest$cluster), ss)
  print(as.vector(round(withinss,4)))
  #tot.withinss <- sum(withinss) # or  
  totss <- ss(x$clusobjbest$obscoord) # or tot.withinss + betweenss
  cat(" (between_SS / total_SS = ",round((betweenss/totss)*100,2),"%)","\n")
  
  cat(paste("\nObjective criterion value:",round(x$clusobjbest$criterion,4),"\n"))
  
  cat("\nAvailable output:\n", 
      sep = "\n")
  print(names(x))
  invisible(x)
  
  
}

ss <- function(x) sum(scale(x, scale = FALSE)^2)
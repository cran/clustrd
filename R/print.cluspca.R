## Define a print method that will be automatically dispatched when print()
## is called on an object of class "cluspca"
print.cluspca <- function(x, ...) {
  
  k = length(x$size)
  if (k == 1)
  {
    d = dim(data.frame(x$attcoord))[2]
    size = x$size
    csize = round((table(x$cluster)/sum(table(x$cluster)))*100,digits=1)
    tt = paste('(',csize,'%)',sep="")
    cs = paste(size, tt, sep = " ", collapse = ", ")
    
    cat(paste("PCA Solution ","in ",d ," dimensions. ","\n", sep = ""))
    
    # cat("\nCluster centroids:\n")
    #  xcent = data.frame(round(x$centroid,4))
    # for (i in 1:k) {
    #    rownames(xcent)[i] = paste("Cluster",i)
    #  }
    #  for (i in 1:ncol(xcent)) {
    #    colnames(xcent)[i] = paste0("Dim.",i)
    #  }
    #  print(xcent)
    attc = data.frame(round(x$attcoord,4))
    cat("\nAttribute scores:\n")
    for (i in 1:ncol(attc)) {
      colnames(attc)[i] = paste0("Dim.",i)
    }
    print(attc)
    
    # cat("\nClustering vector:\n")
    #  print(x$cluster)
    
    
    cat("\nAvailable output:\n", 
        sep = "\n")
    print(names(x))
    invisible(x)
    
  }  else {
    d = dim(x$centroid)[2]
    size = x$size
    csize = round((table(x$cluster)/sum(table(x$cluster)))*100,digits=1)
    #print(csize)
    if (x$center == TRUE) {
      centering = "mean centered" } else {
        centering = "not centered"
      }
    if (x$scale == TRUE) {
      scaling = "standardized" } else {
        scaling = "unstandardized"
      }
    tt = paste('(',csize,'%)',sep="")
    cs = paste(size, tt, sep = " ", collapse = ", ")
    cat(paste("Solution with ",k ," clusters of sizes ", paste(cs, collapse = ", ")," in ",d ," dimensions. ", "Variables were ", centering, " and ", scaling,".", "\n", sep = ""))
    
    cat("\nCluster centroids:\n")
    xcent = data.frame(round(x$centroid,4))
    
    for (i in 1:k) {
      rownames(xcent)[i] = paste("Cluster",i)
    }
    for (i in 1:ncol(xcent)) {
      colnames(xcent)[i] = paste0("Dim.",i)
    }
    print(xcent)
    attc = data.frame(round(x$attcoord,4))
    cat("\nVariable scores:\n")
    for (i in 1:ncol(attc)) {
      colnames(attc)[i] = paste0("Dim.",i)
    }
    print(attc)
    
    cat("\nWithin cluster sum of squares by cluster:\n")
    #resid <- x$obscoord - fitted(x) 
    #tot.withinss <- ss(resid)
    #print(tot.withinss)
    x$centroid = as.matrix(x$centroid)
    betweenss <- ss(x$centroid[x$cluster,]) # or 
    #betweenss <- ss(fitted(x))
    withinss <- sapply(split(as.data.frame(x$obscoord), x$cluster), ss)
    print(as.vector(round(withinss,4)))
    #tot.withinss <- sum(withinss) # or  
    totss <- ss(x$obscoord) # or tot.withinss + betweenss
    cat(" (between_SS / total_SS = ",round((betweenss/totss)*100,2),"%)","\n")
    
    cat(paste("\nObjective criterion value:",round(x$criterion,4),"\n"))
    
    cat("\nAvailable output:\n", 
        sep = "\n")
    print(names(x))
    invisible(x)
    
  }
}
ss <- function(x) sum(scale(x, scale = FALSE)^2)
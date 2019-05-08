tune_cluspcamix <- function(data, nclusrange = 2:7, ndimrange = 2:4, criterion = "asw", dst = "full", alpha = NULL, method = "RKM", center = TRUE, scale = TRUE, rotation = "none", nstart = 100, smartStart = NULL, seed = NULL){
  criterion <- match.arg(criterion, c("asw", "ASW","ch","CH","crit","CRIT"), several.ok = T)[1]
  criterion <- tolower(criterion)
  
  dst <- match.arg(dst, c("full", "FULL","low","LOW","Low","Full"), several.ok = T)[1]
  dst <- tolower(dst)
  
  method <- match.arg(method, c("mixedRKM", "mixedrkm", "mixedrKM", "mixedFKM", 
                                "mixedfkm", "mixedfKM"), several.ok = T)[1]
  method <- toupper(method)

  if (is.null(alpha) == TRUE) {
    if (method == "MIXEDRKM") {
      alpha = 0.5
    }
    else if (method == "MIXEDFKM") 
      alpha = 0
  }
  
  critval = matrix(0,max(length(nclusrange)),max(length(ndimrange)))
  
  m = 1
  n = 1
  for (k in nclusrange) {
    for (d in ndimrange) {
      if (k > d) {
        print(paste('Running for',k,'clusters and',d,'dimensions...'))
        #outcluspca <- cluspca(data=data, nclus = k, ndim = d,alpha = alpha,method = method,  center = center, scale = scale, rotation=rotation, nstart = nstart, smartStart = smartStart, seed = seed)
        outcluspcamix <- cluspcamix(data=data, nclus = k, ndim = d, method = method, alpha = alpha, rotation=rotation, nstart = nstart, smartStart = smartStart, seed = seed)
        
        if (criterion == "asw")
        {
          critval[m,n] <- clusval(outcluspcamix, dst = dst)$asw
        }
        if (criterion == "ch")
        {
          critval[m,n] <- clusval(outcluspcamix, dst = dst)$ch
        }
        
        if (criterion == "crit")
        {
          critval[m,n] <- outcluspcamix$criterion
        }
        
      }
      n = n +1
    }
    n = 1
    m = m +1
  }
  
  #replace 0s with NAs
  critval[critval == 0] <- NA
  if (criterion == "crit")
  {
    indk.best <- which(critval == min(critval,na.rm =TRUE), arr.ind = TRUE)[1]
    indd.best <- which(critval == min(critval,na.rm =TRUE), arr.ind = TRUE)[2]
  } else {
    indk.best <- which(critval == max(critval,na.rm =TRUE), arr.ind = TRUE)[1]
    indd.best <- which(critval == max(critval,na.rm =TRUE), arr.ind = TRUE)[2]
  }
  
  k.best <- nclusrange[indk.best]
  d.best <- ndimrange[indd.best]
  # outcluspcabest = cluspca(data = data, nclus = k.best, ndim = d.best, alpha = alpha, method = method,  center = center, scale = scale, rotation = rotation, nstart = nstart, smartStart = smartStart, seed = seed)
  outcluspcamixbest = cluspcamix(data = data, nclus = k.best, ndim = d.best, method = method, alpha = alpha, rotation = rotation, nstart = nstart, smartStart = smartStart, seed = seed)
  
  rownames(critval) = c(nclusrange)
  colnames(critval) = c(ndimrange)
  
  crit.best = round(critval[indk.best, indd.best],3) 
  crit.grid  = round(critval,3)
  
  crit.grid[is.na(crit.grid)]=''
  crit.grid = as.data.frame(crit.grid)
  
  if (criterion == "asw")
    cluasw = clusval(outcluspcamixbest, dst = dst)$cluasw
  else
    cluasw = NULL
  
  out <- list(clusobjbest = outcluspcamixbest, nclusbest = k.best, ndimbest = d.best, critbest = crit.best, critgrid  = crit.grid, crit = criterion, cluasw = cluasw)
  class(out) = "tuneclus"
  out
}

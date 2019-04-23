tune_clusmca <- function(data, nclusrange = 2:5, ndimrange = 2:4, method = "clusCA", criterion = "asw", dst = "full", alphak = .5, nstart = 100, smartStart = NULL, seed = NULL){
  
  criterion <- match.arg(criterion, c("asw", "ASW","ch","CH","crit","CRIT"), several.ok = T)[1]
  criterion <- tolower(criterion)
  
  dst <- match.arg(dst, c("full", "FULL","low","LOW","Low","Full"), several.ok = T)[1]
  dst <- tolower(dst)
  
  method <- match.arg(method, c("clusCA", "clusca","CLUSCA","CLUSca", "ifcb","iFCB","IFCB","mcak", "MCAk", "MCAK","mcaK"), several.ok = T)[1]
  method <- tolower(method)
  
  if (is.null(alphak) == TRUE)
  { 
    alphak = 0.5
  }
  # outclusmca = list()
  critval = matrix(0,max(length(nclusrange)),max(length(ndimrange)))
  
  m = 1
  n = 1
  for (k in nclusrange) {
    for (d in ndimrange) {
      if (k > d) {
        print(paste('Running for',k,'clusters and',d,'dimensions...'))
        outclusmca <- clusmca(data = data, nclus = k, ndim = d,method = method, alphak = alphak, nstart = nstart,smartStart = smartStart, seed = seed)
        
        if (criterion == "asw")
        {
          critval[m,n] <- clusval(outclusmca, dst = dst)$asw
        }
        if (criterion == "ch")
        {
          critval[m,n] <- clusval(outclusmca, dst = dst)$ch
        }
        
        if (criterion == "crit")
        {
          critval[m,n] <- outclusmca$criterion
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
    if (method != "clusca") {
      indk.best <- which(critval == min(critval,na.rm =TRUE), arr.ind = TRUE)[1]
      indd.best <- which(critval == min(critval,na.rm =TRUE), arr.ind = TRUE)[2]
      #FIX: in case of tie returns the lowest (more parsimonious)
    } else {
      indk.best <- which(critval == max(critval,na.rm =TRUE), arr.ind = TRUE)[1]
      indd.best <- which(critval == max(critval,na.rm =TRUE), arr.ind = TRUE)[2]
      
    }
    
  } else {
    indk.best <- which(critval == max(critval,na.rm =TRUE), arr.ind = TRUE)[1]
    indd.best <- which(critval == max(critval,na.rm =TRUE), arr.ind = TRUE)[2]
  }
  k.best <- nclusrange[indk.best]
  d.best <- ndimrange[indd.best]
  
  outclusmcabest = clusmca(data = data, nclus = k.best, ndim = d.best,method = method, alphak = alphak, nstart = nstart,smartStart = smartStart, seed = seed)
  rownames(critval) = c(nclusrange)
  colnames(critval) = c(ndimrange)
  
  crit.best = round(critval[indk.best, indd.best],3) 
  crit.grid  = round(critval,3)
  
  crit.grid[is.na(crit.grid)]=''
  crit.grid = as.data.frame(crit.grid)
  out <- list(clusobjbest = outclusmcabest, nclusbest = k.best, ndimbest = d.best, critbest = crit.best, critgrid  = crit.grid)
  class(out) = "tuneclus"
  out
}

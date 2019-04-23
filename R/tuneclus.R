tuneclus <- function(data, nclusrange = 3:4, ndimrange = 2:3, method = c("RKM","FKM","mixedRKM","mixedFKM","clusCA","iFCB","MCAk"), criterion = "asw", dst = "full", alpha = NULL, alphak = NULL, center = TRUE, scale = TRUE, rotation = "none", nstart = 100, smartStart = NULL, seed = NULL){
  #wrapper for functions tune_cluspca(), tune_clusmca()
  method <- match.arg(method, c("mixedRKM", "mixedrkm","mixedrKM","mixedFKM", "mixedfkm","mixedfKM","RKM", "rkm","rKM","FKM", "fkm","fKM","clusCA", "clusca","CLUSCA","CLUSca", "ifcb","iFCB","IFCB","mcak", "MCAk", "MCAK","mcaK"), several.ok = T)[1]
  method <- tolower(method)
  
  if (ndimrange[1] >= nclusrange[1]) {
    stop('The number of dimensions must be smaller than the number of clusters.')
  }
  
  if (method %in% c("rkm","fkm")) {
    out = tune_cluspca(data, nclusrange, ndimrange, method = method, criterion = criterion, dst = dst, alpha = alpha, center = center, scale = scale, rotation = rotation, nstart = nstart, smartStart = smartStart, seed = seed)
  } else if (method %in% c("clusca","ifcb","mcak")) {
    out = tune_clusmca(data, nclusrange, ndimrange, method = method, criterion = criterion, dst = dst,  alphak = alphak, nstart = nstart, smartStart = smartStart, seed = seed)
  } else if (method %in% c("mixedrkm","mixedfkm")) {
    out = tune_cluspcamix(data, nclusrange, ndimrange, method = method, criterion = criterion, dst = dst, nstart = nstart, smartStart = smartStart, seed = seed)
  }
  class(out) = "tuneclus"
  out
}
global_bootclus <- function(data, nclusrange = 3:4, ndim = NULL, method = c("RKM","FKM","mixedRKM","mixedFKM","clusCA","MCAk","iFCB"), nboot = 10, alpha = NULL, alphak = NULL, center = TRUE, scale = TRUE, nstart = 100, smartStart = NULL, seed = NULL){
  #wrapper for functions boot_cluspca(), boot_clusmca()
  method <- match.arg(method, c("mixedrkm","MIXEDRKM","mixedRKM","mixedfkm","MIXEDFKM","mixedFKM","RKM", "rkm","rKM","FKM", "fkm","fKM","clusCA", "clusca","CLUSCA","CLUSca", "ifcb","iFCB","IFCB","mcak", "MCAk", "MCAK","mcaK"), several.ok = T)[1]
  method <- tolower(method)
  
  if (!is.null(ndim)) {
    if(ndim >= nclusrange[1]) {
      stop('The number of dimensions must be smaller than the number of clusters.')
    }
  }

  if (method %in% c("rkm","fkm")) {
    out = boot_cluspca(data, krange = nclusrange, nd = ndim, method = method, nboot = nboot,  alpha = alpha, center = center, scale = scale, nstart = nstart, smartStart = smartStart, seed = seed)
  } else if (method %in% c("clusca","ifcb","mcak")) {
    out = boot_clusmca(data, krange = nclusrange, nd = ndim,  method = method, nboot = nboot,  alphak = alphak, nstart = nstart, smartStart = smartStart, seed = seed)
  } else if (method %in% c("mixedrkm","mixedfkm")) {
    out = boot_cluspcamix(data, krange = nclusrange, nd = ndim, nboot = nboot,  alpha = alpha, center = center, scale = scale,  nstart = nstart, smartStart = smartStart, seed = seed)
  }
  class(out) = "genbootclus"
  out
}
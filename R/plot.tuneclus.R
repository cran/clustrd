plot.tuneclus<-function(x, dims = c(1,2), what = c(TRUE,TRUE), cludesc = FALSE, topstdres = 20, attlabs = NULL, binary = FALSE, subplot = FALSE, ...){
  
  out=list()
  if (inherits(x$clusobjbest,"cluspca")) {
    out = plot.cluspca(x$clusobjbest, dims = dims, what = what, cludesc = cludesc, ...)
  }

  if (inherits(x$clusobjbest,"cluspcamix")) {
    out = plot.cluspcamix(x$clusobjbest, dims = dims, what = what, cludesc = cludesc, ...)
  }
  
  if (inherits(x$clusobjbest,"clusmca")) {
    out = plot.clusmca(x$clusobjbest,dims = dims, what = what, cludesc = cludesc, topstdres = topstdres, attlabs = attlabs, binary = binary, subplot = subplot, ...)
  }
  
  out  
}
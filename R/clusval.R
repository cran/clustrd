clusval<-function(obj,dst="full"){
  if(dst=="full"){
    if(class(obj)=="cluspca") {
      data = scale(obj$odata, center = obj$center, scale = obj$scale)
      oDist = daisy(data,metric="euclidean")
    }else{
      oDist=daisy(obj$odata,metric="gower")
    }
  }else{
    oDist=daisy(obj$obscoord,metric="euclidean")
  }
  
  clu_res=cluster.stats(d=oDist,obj$cluster,wgap=F,sepindex=F,sepwithnoise=F)
  
  out=list()
  
  out$ch=clu_res$ch
  out$asw=clu_res$avg.silwidth
  #out$crit=x$criterion
  class(out) = "clusval"
  return(out)
}

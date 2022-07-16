clusval<-function(obj,dst="full"){
  if(dst=="full"){
    if (inherits(obj,"cluspca"))
      {
      data = scale(obj$odata, center = obj$center, scale = obj$scale)
      oDist = daisy(data,metric="euclidean")
    }else{ #clusmca, cluspcamix
      oDist=daisy(obj$odata,metric="gower")
    }
  }else{
    oDist=daisy(obj$obscoord,metric="euclidean")
  }
  
  clu_res=cluster.stats(d=oDist,obj$cluster,wgap=F,sepindex=F,sepwithnoise=F)
  out=list()
  
  out$ch=clu_res$ch
  out$asw=clu_res$avg.silwidth
  out$cluasw=clu_res$clus.avg.silwidths
  #out$crit=x$criterion
  class(out) = "clusval"
  return(out)
}

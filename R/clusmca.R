clusmca <- function(data,nclus,ndim,method=c("clusCA","iFCB","MCAk"),alphak = .5,nstart=100,smartStart=NULL,gamma = TRUE,inboot=FALSE,seed=NULL){

  #### A single cluster gives the MCA solution
  if (nclus == 1) { 
    nstart = 1
    data = data.frame(data, stringsAsFactors = TRUE)
    n = nrow(data)
    #asymmetric map, biplot
    A = mjca(data)$colcoord#[,1:ndim]
    Fm = mjca(data)$rowpcoord#[,1:ndim]
    
    if (gamma == TRUE) {
      distB = sum(diag(t(A)%*%  A))
      g = ((nclus/ncol(data))* distB)^.25 
      A = (1/g)*A
      Fm = g*Fm
    }
    out=list()
    out$obscoord=Fm # observations coordinates 
    out$attcoord=A # attributes coordinates 
    out$centroid = 0#center
    out$cluster = rep(1,n)#cluster
    out$criterion=1 # criterion
    out$size=n #as.integer(aa)  #round((table(cluster)/sum( table(cluster)))*100,digits=1)
    out$odata=data.frame(lapply(data.frame(data),factor),stringsAsFactors = TRUE)
    out$nstart = nstart
    class(out)="clusmca"
    return(out)
  } else {
    
    if (missing(ndim)) {
      warning('The ndim argument is missing. ndim was set to nclus - 1')
      ndim = nclus - 1
    }
    
    if (ndim > nclus) {
      stop('The number of clusters should be larger than the number of dimensions.')
    }
    
    if (ncol(data) < ndim) {
      stop('The number of dimensions should be less than the number of variables.')
    }
    
    method <- match.arg(method, c("clusCA", "clusca","CLUSCA","CLUSca", "ifcb","iFCB","IFCB","mcak", "MCAk", "MCAK","mcaK"), several.ok = TRUE)[1]
    method <- tolower(method)
    
    if(method=="clusca"){
      out=clusCA(data=data,nclus=nclus,ndim=ndim,nstart=nstart,smartStart=smartStart, gamma = gamma,seed=seed,inboot=inboot)
    }
    if(method=="ifcb"){
      out=iFCB(data=data,nclus=nclus,ndim=ndim,nstart=nstart,smartStart=smartStart, gamma = gamma,seed=seed)
    }
    if(method=="mcak"){
      out=MCAk(data=data,nclus=nclus,ndim=ndim,nstart=nstart,alphak = alphak,smartStart=smartStart, gamma = gamma,seed=seed)
    }
    return(out)
  }
}


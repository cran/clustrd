iFCB<- function(data,nclus=3,ndim=2,nstart=100,smartStart=NULL,gamma = TRUE,seed=NULL, inboot = FALSE){
  
  data = data.frame(data)
  minobs = min(sapply(apply(data,2,unique),length))
  maxobs = max(sapply(apply(data,2,unique),length))
  #data=data.frame(data)
  q=ncol(data)
  n=nrow(data)
  
  if (inboot == FALSE) {
    data=as.data.frame(lapply(data,as.factor))
    lab1a=names(data)
    lab1b=lapply(data,function(z) levels(as.factor(z)))
    lab1=rep(lab1a,times=unlist(lapply(lab1b,length)))
    lab2=unlist(lab1b)
    dZ = as.matrix(tab.disjonctif(data))#as.matrix(dummy.data.frame(data, dummy.classes = "ALL")) #as.matrix(dummy.data.frame(data,dummy.classes = "ALL"))
  }
  else 
    dZ = data
  
  ndZ = ncol(dZ)
  
  if(!is.null(seed)) set.seed(seed)
  seed <- round(2^31 * runif(nstart, -1, 1))
  
  best_f=1000000
  pb <- txtProgressBar(style = 3)
  prog = 1
  #  fvec=c()
  for(b in 1:nstart){
    if (b > (prog * (nstart/10))) {
      prog <- prog + 1
    }
    setTxtProgressBar(pb, prog * 0.1)
    # Starting method
    if(is.null(smartStart)){
      #  myseed=seed+b
      #  set.seed(myseed)
      set.seed(seed[b])
      randVec= matrix(ceiling(runif(n)*nclus),n,1)
    }else{
      randVec=smartStart
    }
    
    
    C=tab.disjonctif(randVec)
    
    w= -Inf
    ceps=0.00000001
    itmax=100
    it=0 ### inizializzazione contatore iterazioni
    imp=100000 
    f0 =10000000 ## inizializzazione criterio di arresto
    u=matrix(1,n,1);## vettore di 'uno' 
    #story.obscoord=list()
    
    while((it<=itmax)&&(imp>ceps)){
      
      it=it+1 
      Fmat = t(C) %*% dZ #crossprod(C,dZ) #
      P=Fmat/sum(Fmat)
      r= rowSums(P)
      c= colSums(P) #apply(P,2,sum)
      
      #r=t(t(r))
      #c=t(t(c))
      onec=matrix(1,nrow=ncol(P))
      nsSpc=sqrt(q)*   t(t(t(t(P)*as.vector(1/c)) - r %*% t(onec)) * as.vector(sqrt(c)))
      #
      
      nssvdres=svd(nsSpc)
      
      # print("decomposition done")
      
      nU  = nssvdres$u[,1:ndim]
      nV  = nssvdres$v[,1:ndim]
      nsv = nssvdres$d[1:ndim]
      
      G =  t(t(nU)*nsv)
      
      B =   (1/(sqrt(as.vector(c)*sqrt(q)))) * t(t(nV) * nsv)  
      
      Csize = colSums(C)
      Cw = as.vector(C %*% Csize)
      Y = Cw * (dZ / (n*sqrt(q))) %*% t(t(B) * as.vector(1/nsv))
      
      outK=try(kmeans(Y,centers=G,nstart=100),silent=T)
      if(is.list(outK)==F){
        outK=EmptyKmeans(Y,centers=G)  
        #break
      }
      
      G = outK$centers
      ngvec = outK$cluster
      C = tab.disjonctif(ngvec)
      centerC = matrix(apply(C,2,sum),nrow=n,ncol=nclus,byrow=T)/n
      Zstar = dZ - matrix(c*(n*q)/n,nrow=n,ncol = ndZ,byrow=T)/n
      Cstar = C - centerC
      
      fA=Cstar- t(t(Zstar) * as.vector(sqrt(c))) %*% B %*% t(nU) 
      flossA=sum(diag(t(fA) %*% fA))  
      flossB=sum(diag(t(Y - Cstar %*% G) %*% (Y - Cstar %*% G)))    
      f=flossA+flossB 
      imp=f0-f
      #      fvec=c(fvec,f)
      f0=f
    }
    
    if(f<best_f){
      
      #####gamma scaling
      if (gamma == TRUE) {
        distB = sum(diag(t(B)%*%  B))
        distG = sum(diag(t(G)%*% G))
        g = ((nclus/q)* distB/distG)^.25
        
        B = (1/g)*B
        G = g*G 
        Y = g*Y
      }
      #########################
      
      
      #   best_lam=lambda
      best_f=f
      best_ngvec=ngvec	
      best_Y=Y
      best_B=B
      best_G=G
      best_it=it
    }
    
    
    
  }## end FOR
  
  #####################
  ####################
  #lambda=best_lam
  f=best_f
  ngvec=best_ngvec	
  Y=best_Y	
  B=best_B
  G=best_G
  
  
  it=best_it
  # inert=sum(lambda[1:2]^2)
  #  t_inert=sum(lambda^2)
  
  ####################
  #####################
  if ((minobs==2) & (maxobs==2)) {
    B=B[seq(from=1,to=(2*q),by=2),]
  }
  
  
  #################################  
  #  C=dummy(ngvec)
  #  Dzi=diag(colSums(dZ)^-1)
  #  MZ=scale(dZ,scale=FALSE)
  #  MzDzi= t(C)%*%MZ%*%Dzi
  #  GB=G%*%t(B)
  #  MzDzi_GB=MzDzi-GB
  #  final_loss=sum(diag(t(MzDzi_GB)%*%MzDzi_GB))
  #################################  
  
  cluster = ngvec
  #library(plyr)
  ##reorder cluster membership according to cluster size
  size = table(cluster) #round((table(cluster)/sum( table(cluster)))*100,digits=2)
  aa = sort(size,decreasing = TRUE)
  #aa = sort(csize,decreasing = TRUE)
  cluster = mapvalues(cluster, from = as.integer(names(aa)), to = as.integer(names(table(cluster))))
  #reorder centroids
  G = G[as.integer(names(aa)),]
  setTxtProgressBar(pb, 1)
  out=list() 
  out$obscoord=Y
  rownames(out$obscoord) = rownames(data)
  out$attcoord=B
  rownames(out$attcoord) = paste(lab1,lab2,sep=".")
  out$centroid=G
  cluster = as.integer(cluster)
  names(cluster) = rownames(data) 
  out$cluster=cluster
  #  out$final_loss=final_loss
  out$criterion=f
  #  out$iters=it
  #  out$expl_inertia= (inert/t_inert)
  #  out$lambda=lambda
  # out$critvec=fvec
  out$size=as.integer(aa) #round((table(cluster)/sum( table(cluster)))*100,digits=1)
  out$odata=data.frame(lapply(data.frame(data),factor))
  out$nstart = nstart
  class(out)="clusmca"
  return(out)
}

txtProgressBar <- function(min = 0, max = 1, initial = 0, char = "=", width = NA, 
                           title, label, style = 1, file = "") 
{
  if (!identical(file, "") && !(inherits(file, "connection") && 
                                isOpen(file))) 
    stop("'file' must be \"\" or an open connection object")
  if (!style %in% 1L:3L) 
    style <- 1
  .val <- initial
  .killed <- FALSE
  .nb <- 0L
  .pc <- -1L
  nw <- nchar(char, "w")
  if (is.na(width)) {
    width <- getOption("width")
    if (style == 3L) 
      width <- width - 10L
    width <- trunc(width/nw)
  }
  if (max <= min) 
    stop("must have 'max' > 'min'")
  up1 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    if (.nb < nb) {
      cat(strrep(char, nb - .nb), file = file)
      flush.console()
    }
    else if (.nb > nb) {
      cat("\r", strrep(" ", .nb * nw), "\r", strrep(char, 
                                                    nb), sep = "", file = file)
      flush.console()
    }
    .nb <<- nb
  }
  up2 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    if (.nb <= nb) {
      cat("\r", strrep(char, nb), sep = "", file = file)
      flush.console()
    }
    else {
      cat("\r", strrep(" ", .nb * nw), "\r", strrep(char, 
                                                    nb), sep = "", file = file)
      flush.console()
    }
    .nb <<- nb
  }
  up3 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    pc <- round(100 * (value - min)/(max - min))
    if (nb == .nb && pc == .pc) 
      return()
    cat(paste0("\r  |", strrep(" ", nw * width + 6)), file = file)
    cat(paste(c("\r  |", rep.int(char, nb), rep.int(" ", 
                                                    nw * (width - nb)), sprintf("| %3d%%", pc)), collapse = ""), 
        file = file)
    flush.console()
    .nb <<- nb
    .pc <<- pc
  }
  getVal <- function() .val
  kill <- function() if (!.killed) {
    cat("\n", file = file)
    flush.console()
    .killed <<- TRUE
  }
  up <- switch(style, up1, up2, up3)
  up(initial)
  structure(list(getVal = getVal, up = up, kill = kill), class = "txtProgressBar")
}

setTxtProgressBar <- function (pb, value, title = NULL, label = NULL) 
{
  if (!inherits(pb, "txtProgressBar")) 
    stop(gettextf("'pb' is not from class %s", dQuote("txtProgressBar")), 
         domain = NA)
  oldval <- pb$getVal()
  pb$up(value)
  invisible(oldval)
}


local_bootclus <- function(data, nclus, ndim = NULL, method=c("RKM","FKM","mixedRKM","mixedFKM","clusCA","MCAk","iFCB"), scale = TRUE, center= TRUE, alpha = NULL, nstart=100, nboot=10, alphak = .5, seed=NULL)
{
  clu={}
  if(!is.null(seed)) set.seed(seed)
  seed <- round(2^31 * runif(nboot, -1, 1))
  
  x = data.frame(data,stringsAsFactors = TRUE)
  k = nclus
  nd = ndim
  nk <- length(k)
  nx <- nrow(x)
  
  index1 <- index2 <- indextest <- list()    
  for(b in 1:nboot){
    train1=sample(1:nx, replace=TRUE)
    train2=sample(1:nx, replace=TRUE)
    test=1:nx
    index1[[b]] <- train1
    index2[[b]] <- train2
    indextest[[b]] <- test
  }
  
  valname = as.character(1:k)
  
  method <- match.arg(method, c("mixedRKM", "mixedrkm","mixedrKM","mixedFKM", "mixedfkm","mixedfKM", "RKM", "rkm","rKM","FKM", "fkm","fKM","clusCA", "clusca","CLUSCA","CLUSca", "ifcb","iFCB","IFCB","mcak", "MCAk", "MCAK","mcaK"), several.ok = TRUE)[1]
  method <- toupper(method)
  
  if ((method == "CLUSCA") | (method == "IFCB") | (method == "MCAK")) {
    x = tab.disjonctif(x)#dummy.data.frame(x, dummy.classes = "ALL") # The original super indicator
  }
  
  if (method %in% c("MIXEDRKM","MIXEDFKM")) {
    data = x
   numvars <- sapply(data, is.numeric)
    anynum <- any(numvars)
    catvars <- sapply(data, is.factor)
    anyfact <- any(catvars)
    if (!anynum) 
      cat("\nNo continuous (numeric) variables in data! \n")
    if (!anyfact) 
      cat("\nNo categorical (factor) variables in data! \n")
    ind.sup = NULL
    row.w = NULL #weights of the individuals, could be passed as a parameter, as in FAMD()
    
    if (!anynum) 
      cat("\nNo continuous (numeric) variables in data! Use clusmca() \n")
    
    if (!anyfact) 
      cat("\nNo categorical (factor) variables in data! Use cluspca() \n")
    if (is.null(rownames(data))) 
      rownames(data) = 1:nrow(data)
    if (is.null(colnames(data))) 
      colnames(data) = paste("v", 1:ncol(data), sep = "")
    data <- data.frame(data, stringsAsFactors = TRUE)
    data <- droplevels(data)
    numAct <- which(sapply(data, is.numeric))
    facAct <- which(!sapply(data, is.numeric))
    
    QuantiAct <- as.matrix(data[, numAct, drop = FALSE])
    
#    centre <- moy.ptab(QuantiAct, row.w)
#    QuantiAct <- t(t(QuantiAct) - centre)
#    ecart.type <- ec.tab(QuantiAct, row.w)
#    QuantiAct <- t(t(QuantiAct)/ecart.type)
    QualiAct <-  tab.disjonctif(data[, facAct, drop = FALSE])
    attlabs = c(colnames(QuantiAct),(colnames(QualiAct)))
  #  prop <- colSums(QualiAct * (row.w/sum(row.w)))
  #  QualiAct <- t(t(QualiAct) - moy.ptab(QualiAct, row.w)) #this is centering MZ
  #  QualiAct <- t(t(QualiAct)/sqrt(prop))
    x <- cbind(QuantiAct, QualiAct)  
  }
  
  BFUN <- function(b){
    set.seed(seed[b])
    cat('\n')
    print(paste0("nboot = ", b))
    
    clust1 <- clust2 <- matrix(integer(1), nrow=nx, ncol=nk)
    valmat <- matrix(0, ncol=nk, nrow=length(valname)) 
    for(l in 1:nk)
    {
      # if it's a range of clusters
      if (!is.null(nd))
        ndim = nd
      else
        ndim = k-1
      
      if (method %in% c("CLUSCA","IFCB","MCAK"))  {
        cl1 <- clusmca(x[index1[[b]],,drop=FALSE],nclus=k,ndim=ndim,method = method,nstart=nstart, alphak = alphak, seed = seed, inboot = TRUE)
        cl2 <- clusmca(x[index2[[b]],,drop=FALSE],nclus=k,ndim=ndim,method = method,nstart=nstart, alphak = alphak, seed = seed, inboot = TRUE)
      }
      
      if (method %in% c("RKM","FKM") |  (!is.null(alpha))) {
        cl1 <- cluspca(x[index1[[b]],,drop=FALSE],nclus=k,ndim=ndim,method = method,alpha = alpha,nstart=nstart, scale = scale, center = center, seed = seed)
        cl2 <- cluspca(x[index2[[b]],,drop=FALSE],nclus=k,ndim=ndim,method = method,alpha = alpha,nstart=nstart, scale = scale, center = center, seed = seed)
      }
      
      if (method %in% c("MIXEDRKM","MIXEDFKM")) {
        cl1 <- cluspcamix(x[index1[[b]],,drop=FALSE],nclus=k,ndim=ndim,alpha = alpha,nstart=nstart, scale = scale, center = center, seed = seed, inboot = TRUE)
        cl2 <- cluspcamix(x[index2[[b]],,drop=FALSE],nclus=k,ndim=ndim,alpha = alpha,nstart=nstart, scale = scale, center = center, seed = seed, inboot = TRUE)
        
      }
      
      if (method %in% c("RKM","FKM") | (!is.null(alpha))) {
        clall <- cluspca(x[indextest[[b]],],nclus=k[l],ndim=ndim, method = method, alpha = alpha,nstart=nstart, scale = scale, center = center, seed = seed)
      }
      
      if (method %in% c("CLUSCA","IFCB","MCAK")) 
      {
        clall <- clusmca(x[indextest[[b]],],nclus=k[l],ndim=ndim, method = method, nstart=nstart, alphak = alphak, seed = seed, inboot = TRUE)
      }
      
      if (method %in% c("MIXEDRKM","MIXEDFKM")| (!is.null(alpha)))   {
          clall <- cluspcamix(x[indextest[[b]],],nclus=k[l],ndim=ndim, alpha = alpha,nstart=nstart, scale = scale, center = center,seed = seed, inboot = TRUE)
      }
      
      x1 = data.frame(x[index1[[b]],,drop=FALSE],stringsAsFactors = TRUE)
      gm=apply(x1,2,mean)
      x1$clu = cl1$cluster
      clum=(x1 %>% group_by(clu) %>% summarise_all(funs(mean)))
      
      am = rbind(clum[,-1],gm)
      bm = data.frame(am,stringsAsFactors = TRUE)
      #rownames(bm) = c(paste("C",1:nrow(clum),sep=""),"all")
      
      cl1$centers = as.matrix(bm[1:k,])
      x1 = x1[,-ncol(x1)]
      
      #cl1$centroid
      closest.cluster1 <- function(x) {
        cluster.dist <- apply(cl1$centers, 1, function(y) sqrt(sum((x-y)^2)))
        return(which.min(cluster.dist)[1])
      }
      clust1[,l] <- apply(x, 1, closest.cluster1)
      
      x2 = data.frame(x[index2[[b]],,drop=FALSE],stringsAsFactors = TRUE)
      gm=apply(x2,2,mean)
      
      x2$clu = cl2$cluster
      clum=(x2 %>% group_by(clu) %>% summarise_all(funs(mean)))
      
      am = rbind(clum[,-1],gm)
      bm = data.frame(am,stringsAsFactors = TRUE)
      #rownames(bm) = c(paste("C",1:nrow(clum),sep=""),"all")
      
      cl2$centers = as.matrix(bm[1:k,])
      x2 = x2[,-ncol(x2)]
      
      closest.cluster2 <- function(x) {
        cluster.dist <- apply(cl2$centers, 1, function(y) sqrt(sum((x-y)^2)))
        return(which.min(cluster.dist)[1])
      }
      clust2[,l] <- apply(x, 1, closest.cluster2)
      
      valmat[,l] <- clusterwiseScheme(clall, clust1[,l], clust2[,l])
    }
    list(clust1=clust1, clust2=clust2, valmat=valmat)
  }
  
  z <- lapply(as.list(1:nboot), BFUN)
  #if (multicore)
  #z <- mclapply(as.list(1:nboot), BFUN) 
  
  clust1 <- lapply(z, function(x) x$clust1)
  clust2 <- lapply(z, function(x) x$clust2)
  
  valmat <- unlist(lapply(z, function(x) x$valmat))
  dim(valmat) <- c(length(valname), nk, nboot)  # c(length(valname), nk, nboot)
  dimnames(valmat)[[1]] <- valname#scheme@valname
  dimnames(valmat)[[2]] <- k
  
  #if(verbose) cat("\n")
  
  #  if(dim(x$validation)[2]==1){
  #  } else {
  #    if(is.null(ylab)){
  #      if(is.character(which)) ylab <- which
  #      else ylab <- dimnames(x$validation)[[1]][which]
  #    }
  #    boxplot(as.data.frame(t(x$validation[which,,])), ylab=ylab, ...)
  #  }
  
  out=list()
  out$nclus = k
  out$clust1 = clust1 
  out$clust2 = clust2 
  out$index1 = index1
  out$index2 = index2
#  out$indextest = indextest
  out$Jaccard = t(valmat[,1,])
  
  return(out)
}

clusterwiseScheme <- function(object, c1, c2)
{
  k <-  length(object$size)
  c0 <- object$cluster
  
  z1 <- z2 <- matrix(double(1), nrow=k, ncol=k)
  for(m in 1:k){
    ok0 <- (c0==m)
    for(n in 1:k){
      ok1 <- (c1==n)
      ok2 <- (c2==n)
      z1[m,n] <- sum(ok0&ok1)/sum(ok0|ok1)
      z2[m,n] <- sum(ok0&ok2)/sum(ok0|ok2)
    }
  }
  z1 <- apply(z1, 1, max)
  z2 <- apply(z2, 1, max)
  out = (z1+z2)/2
  return(out)
}

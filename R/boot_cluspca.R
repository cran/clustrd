#Global stability of cluspca() as in Dolnicar & Leisch (2010)
#TODOs: check if it makes sense
# add parallelization as in flexclust()
boot_cluspca <- function(data, krange, nd = NULL, method = "RKM", alpha=NULL,scale = TRUE, center= TRUE,nstart=100, nboot=10,  seed=NULL, ...)
{
  clu={}
  #x = scale(data, center = center, scale = scale)
  x = data.frame(data, stringsAsFactors = TRUE)
  
  if(!is.null(seed)) set.seed(seed)
  seed <- round(2^31 * runif(nboot, -1, 1))
  
  nk <- length(krange)
  nx <- nrow(x)
  
  index1 <- matrix(integer(1), nrow=nx, ncol=nboot)
  index2 <- index1
  ## empirical experiments show parallization does not pay for this
  ## (sample is too fast)
  for(b in 1:nboot){
    index1[,b] <- sample(1:nx, nx, replace=TRUE)
    index2[,b] <- sample(1:nx, nx, replace=TRUE)
  }
  
  BFUN <- function(b){
    set.seed(seed[b])
    cat('\n')
    print(paste0("nboot = ", b))
    clust1 <- clust2 <- matrix(integer(1), nrow=nx, ncol=nk)
    rand <- double(nk)
    conc <- double(nk)
    for(l in 1:nk)
    {
      if(nk>1){
        if (!is.null(nd)) {
          if ((length(nd) >1) & (l==1))  {
            cat('\n')
            print('Warning: the number of dimensions (ndim) must be a single number. Automatically set to the first value in the range.')
          }
          ndim = nd[1]
        }
        else
          ndim = krange[l]-1
        cat('\n')
        print(paste0("Running for ",krange[l]," clusters and ",ndim[1]," dimensions."))
        cl1 <- cluspca(x[index1[,b],,drop=FALSE],nclus=krange[l],ndim=ndim,method = method, nstart=nstart, alpha = alpha,scale = scale, center = center, seed = seed)
        cl2 <- cluspca(x[index2[,b],,drop=FALSE],nclus=krange[l],ndim=ndim,method = method,nstart=nstart, alpha = alpha, scale = scale, center = center, seed = seed)
      } else{
        if (!is.null(nd)) {
          if ((length(nd) >1) & (l==1))  {
            cat('\n')
            print('Warning: the number of dimensions (ndim) must be a single number. Automatically set to the first value in the range.')
          }
          ndim = nd[1]
        }
        else
          ndim = krange-1
        cat('\n')
        print(paste0("Running for ",krange," clusters and ",ndim[1]," dimensions."))
        cl1 <- cluspca(x[index1[,b],,drop=FALSE],nclus=krange,ndim=ndim,method = method,nstart=nstart, alpha = alpha, scale = scale, center = center, seed = seed)
        cl2 <- cluspca(x[index2[,b],,drop=FALSE],nclus=krange,ndim=ndim,method = method,nstart=nstart, alpha = alpha,scale = scale, center = center, seed = seed)
      }
      # clall <- cluspca(x,nclus=krange[l],ndim=ndim, method = method, nstart=nstart, alpha = alpha, scale = scale, center = center, seed = seed)
      
      x1 = x[index1[,b],,drop=FALSE]
      gm=apply(x1,2,mean)
      x1$clu = cl1$cluster
      clum=(x1 %>% group_by(clu) %>% summarise_all(list(mean)))
      
      am = rbind(clum[,-1],gm)
      bm = data.frame(am, stringsAsFactors = TRUE)
      #rownames(bm) = c(paste("C",1:nrow(clum),sep=""),"all")
      
      cl1$centers = as.matrix(bm[1:krange[l],])
      x1 = x1[,-ncol(x1)]
      
      #cl1$centroid
      closest.cluster1 <- function(x) {
        cluster.dist <- apply(cl1$centers, 1, function(y) sqrt(sum((x-y)^2)))
        return(which.min(cluster.dist)[1])
      }
      clust1[,l] <- apply(x, 1, closest.cluster1)
      # clall$obscoord
      # print(x)
      x2 = x[index2[,b],,drop=FALSE]
      gm=apply(x2,2,mean)
      
      x2$clu = cl2$cluster
      clum=(x2 %>% group_by(clu) %>% summarise_all(list(mean)))
      
      am = rbind(clum[,-1],gm)
      bm = data.frame(am, stringsAsFactors = TRUE)
      #rownames(bm) = c(paste("C",1:nrow(clum),sep=""),"all")
      
      cl2$centers = as.matrix(bm[1:krange[l],])
      x2 = x2[,-ncol(x2)]
      #cl2$centroid
      closest.cluster2 <- function(x) {
        cluster.dist <- apply(cl2$centers, 1, function(y) sqrt(sum((x-y)^2)))
        return(which.min(cluster.dist)[1])
      }
      #clall$obscoord
      clust2[,l] <- apply(x, 1, closest.cluster2)
      #replace this
      
      rand[l] <- randIndex(clust1[,l], clust2[,l])
      
      I = length(unique(clust1[,l]))
      J = length(unique(clust2[,l]))
      chisq <- suppressWarnings(chisq.test(table(clust1[,l],clust2[,l]))$statistic)
      conc[l]<- chisq/(nx*(sqrt(I*J)-1))
      
      
    }
    list(clust1=clust1, clust2=clust2, rand=rand, conc=conc)
    
  }
  
  ## empirical experiments show parallization does not pay for the 
  ## following (element extraction from list is too fast)
  # print(system.time({
  # z <- mclapply(as.list(1:nboot), BFUN)
  #  }))
  
  #  print(system.time({
  z <- lapply(as.list(1:nboot), BFUN)
  #  }))
  
  clust1 <- unlist(lapply(z, function(x) x$clust1))
  clust2 <- unlist(lapply(z, function(x) x$clust2))
  dim(clust1) <- dim(clust2) <- c(nx, nk, nboot)
  
  if(nk > 1) {
    rand <- t(sapply(z, function(x) x$rand))
    conc <- t(sapply(z, function(x) x$conc))
  }
  else {
    rand <- as.matrix(sapply(z, function(x) x$rand))
    conc <- as.matrix(sapply(z, function(x) x$conc))
  }
  colnames(rand) <- krange
  colnames(conc) <- krange
  
  
  out=list()
  out$nclusrange = krange
  out$clust1 = clust1 
  out$clust2 = clust2 
  out$index1 = index1
  out$index2 = index2
  out$rand = rand
  out$moc = conc
  
  return(out)
  
}

randIndex <- function (x, y) 
{
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y)) 
    stop("arguments must be vectors of the same length")
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1))) 
    return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + 
                                                     a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}

# setGeneric("randIndex", function(x, y, correct=TRUE, original=!correct)
#   standardGeneric("randIndex"))
# 
# setMethod("randIndex", signature(x="ANY", y="ANY"),
#           function(x, y, correct=TRUE, original=!correct){
#             if(correct)
#               comPart(x, y, type="ARI")
#             else
#               comPart(x, y, type="RI")
#           })
# 
# setMethod("randIndex", signature(x="table", y="missing"),
#           doRandIndex <- function(x, y, correct=TRUE, original=!correct)
#           {
#             if(length(dim(x))!=2)
#               stop("Argument x needs to be a 2-dimensional table.")
#             
#             n <- sum(x)
#             ni <- apply(x, 1, sum)
#             nj <- apply(x, 2, sum)
#             n2 <- choose(n, 2)
#             
#             rand <- NULL
#             if(correct){
#               nis2 <- sum(choose(ni[ni > 1], 2))
#               njs2 <- sum(choose(nj[nj > 1], 2))
#               rand <- c(ARI=c(sum(choose(x[x > 1], 2)) -
#                                 (nis2 * njs2)/n2)/((nis2 + njs2)/2 - (nis2 * njs2)/n2))
#             }
#             
#             if(original){
#               rand <- c(rand, RI = 1 + (sum(x^2) - (sum(ni^2) + sum(nj^2))/2)/n2)
#             }
#             
#             return(rand)
#           })
# 
# ###**********************************************************
# 
# setGeneric("comPart", function(x, y, type=c("ARI","RI","J","FM"))
#   standardGeneric("comPart"))
# 
# setMethod("comPart", signature(x="flexclust", y="flexclust"),
#           function(x, y, type){
#             doComPart(clusters(x), clusters(y), type)
#           })
# 
# setMethod("comPart", signature(x="flexclust", y="numeric"),
#           function(x, y, type){
#             doComPart(clusters(x), y, type)
#           })
# 
# setMethod("comPart", signature(x="numeric", y="flexclust"),
#           function(x, y, type){
#             doComPart(x, clusters(y), type)
#           })
# 
# setMethod("comPart", signature(x="numeric", y="numeric"),
#           doComPart <- function(x, y, type=c("ARI","RI","J","FM"))
#           {
#             type <- toupper(type)
#             if(length(x)!=length(y))
#               stop("x an y must have the same length")
#             
#             nxx <- countPairs(x, y)
#             
#             res <- NULL
#             if("ARI" %in% type)
#               res <- c(doRandIndex(table(x,y), correct=TRUE))
#             
#             if("RI" %in% type)
#               res <- c(res, RI=sum(diag(nxx))/sum(nxx))
#             
#             if("J" %in% type)
#               res <- c(res, J=nxx[2,2]/sum(nxx[-1]))
#             
#             if("FM" %in% type){
#               tab <- table(x)
#               w <- sum(tab*(tab-1))/2
#               tab <- table(y)
#               w <- w*sum(tab*(tab-1))/2
#               res <- c(res, FM=nxx[2,2]/sqrt(w))
#             }
#             res        
#           })
# 
# 
# countPairs <- function(x, y)
# {
#   if(length(x)!=length(y))
#     stop("x an y must have the same length")
#   
#   res <- .C(C_countPairs,
#             as.integer(x),
#             as.integer(y),
#             as.integer(length(x)),
#             res=double(4))[["res"]]
#   matrix(res, nrow=2, dimnames=list(0:1,0:1))
# }
# 

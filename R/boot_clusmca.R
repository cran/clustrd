boot_clusmca <- function(data, krange, nd=NULL, method = "clusCA", nstart=100, nboot=10, seed=NULL,...)
{
  clu={}
  #bootstrapping on Z
  data = data.frame(data)
  data=as.data.frame(lapply(data,as.factor))
  x = data.frame(tab.disjonctif(data))
  
  #dummy.data.frame(data, dummy.classes = "ALL")
  if(!is.null(seed)) set.seed(seed)
  seed <- round(2^31 * runif(nboot, -1, 1))
  
  nk <- length(krange)
  nx <- nrow(x)
  
  index1 <- matrix(integer(1), nrow=nx, ncol=nboot)
  index2 <- index1
  
  for(b in 1:nboot){
    index1[,b] <- sample(1:nx, nx, replace=TRUE)
    index2[,b] <- sample(1:nx, nx, replace=TRUE)
  }
  
  BFUN <- function(b){
    
    set.seed(seed[b])
    cat('\n')
    print(paste0("nboot = ", b))
    
    clust1 <- clust2 <- matrix(integer(1), nrow=nx, ncol=nk)
    cent1 <- cent2 <- list()
    rand <- double(nk)
    conc <- double(nk)
    for(l in 1:nk)
    {
      if(nk>1){
        if (!is.null(nd)) {
          if ((length(nd) >1) & (l==1))  {
            cat('\n')
            print('Warning: the number of dimensions (nd) must be a single number. Automatically set to the first value in the range.')
          }
        ndim = nd[1]
      }
      else
        ndim = krange[l]-1
      cat('\n')
      print(paste0("Running for ",krange[l]," clusters and ",ndim[1]," dimensions."))
      
      x1 = x[index1[,b],,drop=FALSE]
      x2 = x[index2[,b],,drop=FALSE]
      
      cl1 <- clusmca(x[index1[,b],,drop=FALSE],nclus=krange[l],ndim=ndim,method = method,nstart=nstart, seed = seed)
      cl2 <- clusmca(x[index2[,b],,drop=FALSE],nclus=krange[l],ndim=ndim,method = method,nstart=nstart, seed = seed)
      
    } else{
      if (!is.null(nd)) {
        if ((length(nd) >1) & (l==1))  {
          cat('\n')
          print('Warning: the number of dimensions must be a single number, not a range. Automatically set to the first value.')
        }
        ndim = nd[1]
      }
      else{
        ndim = krange-1
      }
      cat('\n')
      print(paste0("Running for ",krange," clusters and ",ndim[1]," dimensions."))
      x1 = x[index1[,b],,drop=FALSE]
      x2 = x[index2[,b],,drop=FALSE]
      
      cl1 <- clusmca(x1,nclus=krange,ndim=ndim,method = method,nstart=nstart, seed = seed, inboot = TRUE)
      cl2 <- clusmca(x2,nclus=krange,ndim=ndim,method = method,nstart=nstart, seed = seed, inboot = TRUE)
    }
    
    gm=apply(x1,2,mean)
    x1$clu = cl1$cluster
    clum=(x1 %>% group_by(clu) %>% summarise_all(list(mean)))
    bm = data.frame(rbind(clum[,-1],gm))
    #rownames(bm) = c(paste("C",1:nrow(clum),sep=""),"all")
    cl1$centers = as.matrix(bm[1:krange[l],])
    x1 = x1[,-ncol(x1)]
    
    closest.cluster1 <- function(x) {
      cluster.dist <- apply(data.frame(cl1$centers), 1, function(y) sqrt(sum((x-y)^2)))
      return(which.min(cluster.dist)[1])
    }
    clust1[,l] <- apply(x, 1, closest.cluster1)
    
    #   x2 = dummy.data.frame(x1, dummy.classes = "ALL") # The original super indicator
    gm=apply(x2,2,mean)
    x2$clu = cl2$cluster
    clum=(x2 %>% group_by(clu) %>% summarise_all(list(mean)))
    bm = data.frame(rbind(clum[,-1],gm))
    cl2$centers = as.matrix(bm[1:krange[l],])
    x2 = x2[,-ncol(x2)]
    
    closest.cluster2 <- function(x) {
      cluster.dist <- apply(data.frame(cl2$centers), 1, function(y) sqrt(sum((x-y)^2)))
      return(which.min(cluster.dist)[1])
    }
    clust2[,l] <- apply(x, 1, closest.cluster2)
    
    #     if (measure == "ari") {
    rand[l] <- randIndex(clust1[,l], clust2[,l])
    #      }
    #     if (measure =="conc") {
    I = length(unique(clust1[,l]))
    J = length(unique(clust2[,l]))
    chisq = suppressWarnings(chisq.test(table(clust1[,l],clust2[,l]))$statistic)
    conc[l]<- chisq/(nx*(sqrt(I*J)-1))
    
    #      }
    
    #   if(nrow(cl1@centers) < k[l]) {
    #      extra <- matrix(NA, 
    #                     ncol=ncol(cl1@centers), 
    #                      nrow=k[l]-nrow(cl1@centers))
    #     cent1[[l]] <- rbind(cl1@centers, extra)
    #    }
    
    #    if(nrow(cl2@centers) < k[l]) {
    #      extra <- matrix(NA, 
    #                      ncol=ncol(cl2@centers), 
    #                     nrow=k[l]-nrow(cl2@centers))
    #      cent2[[l]] <- rbind(cl2@centers, extra)
    #    }
  }
  list(clust1=clust1, clust2=clust2, rand=rand,conc=conc)
  
  #  list(cent1=cent1, cent2=cent2, clust1=clust1, clust2=clust2,
  #       rand=rand)
  
}

## empirical experiments show parallization does not pay for the 
## following (element extraction from list is too fast)
#z <- MClapply(as.list(1:nboot), BFUN, multicore=multicore)

z <- lapply(as.list(1:nboot), BFUN)

clust1 <- unlist(lapply(z, function(x) x$clust1))
clust2 <- unlist(lapply(z, function(x) x$clust2))
dim(clust1) <- dim(clust2) <- c(nx, nk, nboot)

#  cent1 <- cent2 <- list()
#  for(l in 1:nk){
#    cent1[[l]] <- unlist(lapply(z, function(x) x$cent1[[l]]))
#    cent2[[l]] <- unlist(lapply(z, function(x) x$cent2[[l]]))
#    dim(cent1[[l]]) <- dim(cent2[[l]]) <- c(k[l], ncol(x), nboot)
#  }

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


# new("bootFlexclust", k=as.integer(k), centers1=cent1, centers2=cent2,
#      cluster1=clust1, cluster2=clust2, index1=index1, index2=index2,
#      rand=rand, call=MYCALL)
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


# setMethod("show", "bootFlexclust",
#           function(object){
#             cat("An object of class", sQuote(class(object)),"\n\n")
#             cat("Call:\n")
#             print(object@call)
#             cat("\nNumber of bootstrap pairs:", nrow(object@rand),"\n")
#           })
# 
# setMethod("summary", "bootFlexclust",
#           function(object){
#             cat("Call:\n")
#             print(object@call)
#             cat("\nSummary of Rand Indices:\n")
#             print(summary(object@rand))
#           })
# 
# setMethod("plot", signature("bootFlexclust","missing"),
#           function(x, y, ...){
#             boxplot(x, ...)
#           })
# 
# setMethod("boxplot", "bootFlexclust",
#           function(x, ...){
#             boxplot(as.data.frame(x@rand), ...)
#           })
# 
# setMethod("densityplot", "bootFlexclust",
#           function(x, data, ...){
#             Rand <- as.vector(x@rand)
#             k <- rep(colnames(x@rand), rep(nrow(x@rand), ncol(x@rand)))
#             k <- factor(k, levels=colnames(x@rand))
#             
#             densityplot(~Rand|k, as.table=TRUE, to=1, ...)
#           })



disp <- function(x, clus, square = TRUE) {
  n <- length(clus)
  k <- max(clus)
  clus <- as.numeric(clus)
  x <- as.matrix(x)
  centers <- matrix(nrow = k, ncol = ncol(x))
  for (i in 1:k) {
    tryCatch(centers[i, ] <- apply(x[clus == i, ], 2, mean), 
             error = function(e) {print(dim(x))})
  }
  sumsq <- rep(0, k)
  if(square == TRUE) 
    x <- (x - centers[clus, ])^2
  else
    x <- abs((x - centers[clus, ]))
  for (i in 1:k) {
    sumsq[i] <- sum(x[clus == i, ])
  }
  sumsq
}

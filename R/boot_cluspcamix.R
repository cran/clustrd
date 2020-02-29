#Global stability of cluspcamix() as in Dolnicar & Leisch (2010)
#TODOs: 
# add parallelization as in flexclust()
boot_cluspcamix <- function(data, krange, nd = NULL,alpha=NULL, scale = TRUE, center= TRUE,nstart=100, nboot=10, seed=NULL, ...)
{
  clu={}
  data = data.frame(data, stringsAsFactors = TRUE)
  numvars <- sapply(data, is.numeric)
  anynum <- any(numvars)
  catvars <- sapply(data, is.factor)
  anyfact <- any(catvars)
  if (!anynum) 
    cat("\nNo continuous (numeric) variables in data! \n")
  if (!anyfact) 
    cat("\nNo categorical (factor) variables in data! \n")
  if (is.null(rownames(data))) 
    rownames(data) = 1:nrow(data)
  if (is.null(colnames(data))) 
    colnames(data) = paste("V", 1:ncol(data), sep = "")
  data <- as.data.frame(data)
  data <- droplevels(data)
  numAct <- which(sapply(data, is.numeric))
  facAct <- which(!sapply(data, is.numeric))
  
  QuantiAct <- as.matrix(data[, numAct, drop = FALSE])
  #numobs = nrow(data)
  #standardize continuous
  # QuantiAct <- t(t(QuantiAct) - as.vector(crossprod(rep(1,numobs)/numobs, as.matrix(QuantiAct))))
  # QuantiAct <- t(t(QuantiAct)/sqrt(as.vector(crossprod(rep(1,numobs)/numobs, 
  #    as.matrix(QuantiAct^2)))))
  
  QualiAct <-  tab.disjonctif(data[, facAct, drop = FALSE])
  #attlabs = c(colnames(QuantiAct),(colnames(QualiAct)))
  
  #standardize categorical
  #  prop <- colSums(QualiAct * (rep(1,numobs)/numobs))
  #  QualiAct <- t(t(QualiAct) - as.vector(crossprod(rep(1,numobs)/numobs, as.matrix(QualiAct)))  ) #this is centering MZ
  #  QualiAct <- t(t(QualiAct)/sqrt(prop))
  QualiAct = data.frame(QualiAct, stringsAsFactors = TRUE)
  for (i in 1:ncol(QualiAct)) 
    QualiAct[,i] = factor((QualiAct[,i]))
  
  x <- cbind(QuantiAct, QualiAct) 
  xgood <- x
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
            print('Warning: the number of dimensions (nd) must be a single number. Automatically set to the first value in the range.')
          }
          ndim = nd[1]
        }
        else
          ndim = krange[l]-1
        cat('\n')
        print(paste0("Running for ",krange[l]," clusters and ",ndim[1]," dimensions."))
        cl1 <- cluspcamix(x[index1[,b],,drop=FALSE],nclus=krange[l],ndim=ndim,nstart=nstart, scale = scale, center = center, seed = seed, inboot = TRUE)
        cl2 <- cluspcamix(x[index2[,b],,drop=FALSE],nclus=krange[l],ndim=ndim,nstart=nstart, scale = scale, center = center, seed = seed, inboot = TRUE)
      } else{
        if (!is.null(nd)) {
          if ((length(nd) >1) & (l==1))  {
            cat('\n')
            print('Warning: the number of dimensions (nd) must be a single number. Automatically set to the first value in the range.')
          }
          ndim = nd[1]
        }
        else
          ndim = krange-1
        cat('\n')
        print(paste0("Running for ",krange," clusters and ",ndim[1]," dimensions."))
        cl1 <- cluspcamix(x[index1[,b],,drop=FALSE],nclus=krange,ndim=ndim,nstart=nstart, scale = scale, center = center, seed = seed, inboot = TRUE)
        cl2 <- cluspcamix(x[index2[,b],,drop=FALSE],nclus=krange,ndim=ndim,nstart=nstart, scale = scale, center = center, seed = seed, inboot = TRUE)
      }
      
      indx <- sapply(xgood, is.factor)
      xgood[indx] <- lapply(x[indx], function(x) as.numeric(as.character(x)))
      x1 = xgood[index1[,b],]
      gm=apply(x1,2,mean)
      x1$clu = cl1$cluster
      clum=(x1 %>% group_by(clu) %>% summarise_all(list(mean)))
      
      bm = data.frame(rbind(clum[,-1],gm),stringsAsFactors = TRUE)
      #rownames(bm) = c(paste("C",1:nrow(clum),sep=""),"all")
      cl1$centers = as.matrix(bm[1:krange[l],])
      x1 = x1[,-ncol(x1)]
      closest.cluster1 <- function(x) {
        cluster.dist <- apply(data.frame(cl1$centers, stringsAsFactors = TRUE), 1, function(y) sqrt(sum((x-y)^2)))
        return(which.min(cluster.dist)[1])
      }
      clust1[,l] <- apply(xgood, 1, closest.cluster1)
      
      
      x2 = xgood[index2[,b],,drop=FALSE]
      gm=apply(x2,2,mean)
      
      x2$clu = cl2$cluster
      clum=(x2 %>% group_by(clu) %>% summarise_all(list(mean)))
      
      bm = data.frame(rbind(clum[,-1],gm), stringsAsFactors = TRUE)
      
      cl2$centers = as.matrix(bm[1:krange[l],])
      x2 = x2[,-ncol(x2)]
      
      closest.cluster2 <- function(x) {
        cluster.dist <- apply(data.frame(cl2$centers, stringsAsFactors = TRUE), 1, function(y) sqrt(sum((x-y)^2)))
        return(which.min(cluster.dist)[1])
      }
      clust2[,l] <- apply(xgood, 1, closest.cluster2)
      
      rand[l] <-  randIndex(clust1[,l], clust2[,l])
      
      I = length(unique(clust1[,l]))
      J = length(unique(clust2[,l]))
      chisq = suppressWarnings(chisq.test(table(clust1[,l],clust2[,l]))$statistic)
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




###**********************************************************

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

MCAk <- function(data, nclus = 3, ndim = 2, alphak = .5, nstart = 100, smartStart=NULL,gamma = TRUE, seed=NULL){
  out=list()
  
  group={}
  data=data.frame(data)
  if(!is.null(seed)) set.seed(seed)
  seed <- round(2^31 * runif(nstart, -1, 1))
  
  if (alphak == 1) { #Tandem approach MCA + k-means
    
    n = nrow(data)
    #asymmetric map, biplot
    A = mjca(data)$colcoord[,1:ndim]
    Fm = mjca(data)$rowpcoord[,1:ndim]
    
    if(is.null(smartStart)){
      set.seed(seed[1])
      randVec= matrix(ceiling(runif(n)*nclus),n,1)
    }else{
      randVec=smartStart
    }
    
    U = dummy(randVec)
    
    center = chol2inv(chol(crossprod(U))) %*% t(U) %*% Fm 
    outK = try(kmeans(Fm, centers = center, nstart = 100), silent = TRUE)
    
    center=outK$centers
    
    if(is.list(outK) == FALSE){
      outK = EmptyKmeans(Fm,centers=center)  
      # break
    }
    
    index = outK$cluster
    cluster = as.numeric(index)
    cluster = as.integer(cluster)
    names(cluster) = rownames(data) 
    size = table(cluster)
    aa = sort(size,decreasing = TRUE)
    
    if (gamma == TRUE) {
      distB = sum(diag(crossprod(A)))
      
      distG = sum(diag(crossprod(center)))
      g = ((nclus/ncol(data))* distB/distG)^.25
      
      A = (1/g)*A
      center = g*center #is this needed
      Fm = g*Fm
    }
    
    cluster = mapvalues(cluster, from = as.integer(names(aa)), to = as.integer(names(table(cluster))))
    #reorder centroids
    center = center[as.integer(names(aa)),]
    out$obscoord=Fm # observations coordinates 
    out$attcoord=A # attributes coordinates 
    out$centroid = center
    out$cluster = cluster
    out$criterion = 1 # criterion
    out$size=as.integer(aa)  #round((table(cluster)/sum( table(cluster)))*100,digits=1)
    out$odata=data.frame(lapply(data.frame(data),factor))
    out$nstart = nstart
    class(out)="clusmca"
    return(out)
  } else {
    
    zz = as.matrix(dummy.data.frame(data,dummy.classes = "ALL"))
    # data = data.matrix(data)
    data=data.frame(data)
    n = nrow(data)
    zitem = ncol(data)            
    zncati = sapply(data.frame(data), function(x) length(unique(x))) #apply(data,2,max)
    oner = matrix(1,n,1)
    muz  = colMeans(zz)
    z = zz - oner %*% muz  
    rr = colSums(crossprod(zz))
    invsqDr = 1/sqrt(rr) 
    #  invsqDr =diag(1/sqrt(rr))
    
    M = z
    Pz = t(t(z) * as.vector(invsqDr) )  
    PPz = crossprod(Pz)
    
    svdPz = svd(PPz)
    invsqD = 1/sqrt(svdPz$d[1:ndim])
    F0 = t(t(Pz %*%  svdPz$v[,1:ndim]) * as.vector(invsqD))
    
    oldf = 1e+06
    pb <- txtProgressBar(style = 3)
    prog = 1
    
    for(b in 1:nstart){
      if (b > (prog * (nstart/10))) {
        prog <- prog + 1
      }
      setTxtProgressBar(pb, prog * 0.1)
      
      Fv={}
      Fm = F0
      # Starting method
      if(is.null(smartStart)){
       # myseed = seed+b
      #  set.seed(myseed)
        set.seed(seed[b])
        randVec= matrix(ceiling(runif(n)*nclus),n,1)
      }else{
        randVec = smartStart
      }
      
      U = dummy(randVec)
      
      #   center = pseudoinverse(t(U) %*% U) %*% t(U) %*% Fm 
      
      mydata = as_tibble(cbind(Fm,group = as.factor(randVec)))
      
      center=mydata%>%
        group_by(group) %>%
        summarise_all(mean) #%>%
    #    select(-group)
      
      center = center[,-1]
      itmax = 100
      it = 0
      ceps = 1e-04
      imp = 1e+05
      f0 = 1e+05
      
      while((it <= itmax ) && ( imp > ceps ) ){
        it=it+1
        
        ####################################################################
        ## STEP 1: update of U #############################################
        ####################################################################
        #use Lloyd's k-means algorithm to get the results of Hwang and Takane (2006)
        #outK=try(kmeans(Fm,centers=center,algorithm="Lloyd",nstart=100),silent=T)
        outK = try(kmeans(Fm,centers=center,nstart=100),silent=T)
        
        if(is.list(outK)==F){
          outK = EmptyKmeans(Fm,centers=center)  
          #  break
        }
        center=outK$centers
        index = outK$cluster
        U = dummy(index)
        U0 = scale(U,center=TRUE, scale=FALSE)
        uu = colSums(crossprod(U)) #colSums(t(U)%*% U)
        invsqDru = 1/sqrt(c(rr,uu))
        
        ####################################################################
        ## STEP 2: update of Fm and Wj ######################################
        ####################################################################
        MU = cbind(alphak*M,(1-alphak)*U0)
        Pzu = t(t(MU) * as.vector(invsqDru))
        
        Pzu[is.nan(Pzu)] <- 0
        PPzu = t(Pzu) %*% Pzu
        svdPzu = svd(PPzu)
        #invsqD = diag(1/sqrt(svdPzu$d))
        
        Fm = t(t(Pzu %*% svdPzu$u[,1:ndim]) * as.vector(1/sqrt(svdPzu$d[1:ndim])))
        ft1 = 0
        k = 1
        kk = 0
        kk = kk+zncati[1]
        Tm = z[,k:kk]
        #chol2inv(chol(crossprod(Tm))) gives an error!
        W = pseudoinverse(crossprod(Tm))%*% t(Tm) %*% Fm 
        A = W
        ft1 = ft1+sum(diag((crossprod(Fm))-(t(Fm) %*% Tm %*% W)))
        k = kk+1
        for(j in 2:zitem){
          kk = kk+zncati[j]
          Tm = z[,k:kk]
          W = pseudoinverse(crossprod(Tm))%*% t(Tm) %*% Fm 
          A = rbind(A,W)
          ft1 = ft1 + sum(diag((crossprod(Fm))-(t(Fm) %*% Tm %*% W))) ## MCA
          k=kk+1
        }
        ft2 = sum(diag((crossprod(Fm)) - (t(Fm) %*% U %*% center)))
        #check again
        f =  alphak*ft1 + (1-alphak)*ft2
        
        imp=f0-f
        f0=f
        Fv = cbind(Fv,f)
      } # end WHILE
      
      if (f < oldf){
        #####gamma scaling
        if (gamma == TRUE) {
          distB = sum(diag(crossprod(A)))
          distG = sum(diag(crossprod(center)))
          g = ((nclus/zitem)* distB/distG)^.25
          
          A = (1/g)*A
          center = g*center #is this needed
          Fm = g*Fm
        }
        #########################
        oldF = Fm
        oldindex = index
        oldf = f				
        Uold = U				
        Aold = A
        centerold = center
      }
    } ##end of FOR
    
    
    Fm = oldF
    index = oldindex
    cluster = as.numeric(index)
    f = oldf
    U = Uold
    A = Aold
    # it=itold
    center = centerold
    
    size = table(cluster)
    aa = sort(size,decreasing = TRUE)
    
    cluster = mapvalues(cluster, from = as.integer(names(aa)), to = as.integer(names(table(cluster))))
    #reorder centroids
    center = center[as.integer(names(aa)),]
    setTxtProgressBar(pb, 1)
    
    out$obscoord = Fm # observations coordinates 
    out$attcoord = A # attributes coordinates 
    out$centroid = center # centroids
    cluster = as.integer(cluster)
    names(cluster) = rownames(data) 
    out$cluster = cluster #as.numeric(index) # cluster membership
    out$criterion = f # criterion
    out$size = as.integer(aa)  #round((table(cluster)/sum( table(cluster)))*100,digits=1)
    out$odata = data.frame(lapply(data.frame(data),factor))
    out$nstart = nstart
    class(out) = "clusmca"
    return(out)
  }  
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


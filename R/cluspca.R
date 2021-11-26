cluspca <- function(data, nclus, ndim, alpha=NULL, method=c("RKM","FKM"), center = TRUE, scale = TRUE, rotation="none", nstart=100, smartStart=NULL, seed=NULL)
{
  # Code optimized for speed by 2KC - June 2019
  
  group=trueOrd=clu={}
  . = NULL
  if (nrow(data) < 700) { #threshold for largish data sets
    #### A single cluster gives the PCA solution
    if (nclus == 1) {
      nstart = 1
      data = data.frame(data, stringsAsFactors = TRUE)
      n = nrow(data)
      #asymmetric map, biplot
      outp = princomp(data)
      out=list()
      out$obscoord=outp$scores[,1:ndim] # observations coordinates
      out$attcoord=data.matrix(outp$loadings[,1:ndim]) # attributes coordinates
      rownames(out$obscoord) = rownames(data)
      rownames(out$attcoord) = colnames(data)
      
      out$centroid = 0 #center
      out$cluster = rep(1,n)#cluster
      names(out$cluster) = rownames(data)
      out$criterion = 1 # criterion
      out$size=n #as.integer(aa)  #round((table(cluster)/sum( table(cluster)))*100,digits=1)
      out$odata= data.frame(data, stringsAsFactors = TRUE)#data.frame(lapply(data.frame(data),factor))
      out$scale = scale
      out$center = center
      out$nstart = nstart
      class(out)="cluspca"
      return(out)
    } else {
      
      #NOTE: FactorialKM needs smartstart k-means or else to perform well
      if (missing(ndim)) {
        warning('The ndim argument is missing. ndim was set to nclus - 1')
        ndim = nclus - 1
      }
      
      if (ndim > nclus) {
        stop('The number of clusters should be more than the number of dimensions.')
      }
      
      if (ncol(data) < ndim) {
        stop('The number of dimensions should be less than the number of variables.')
      }
      
      method <- match.arg(method, c("RKM", "rkm","rKM","FKM", "fkm","fKM"), several.ok = TRUE)[1]
      method <- toupper(method)
      
      if(!is.null(seed)) set.seed(seed)
      seed <- round(2^31 * runif(nstart, -1, 1))
      #  If alpha = .5 gives RKM, alpha=1 PCA and alpha =0  FKM.
      if (is.null(alpha) == TRUE)
      {
        if (method == "RKM") {
          alpha = .5
        } else if (method == "FKM") {
          alpha = 0
        }
      }
      odata = data
      data =  scale(data, center = center, scale = scale)
      
      data = data.matrix(data)
      n = dim(data)[1]
      m = dim(data)[2]
      conv=1e-6  # convergence criterion
      bestf = 10^12
      pb <- txtProgressBar(style = 3)
      prog = 1
      for (run in c(1:nstart)) {
        if (run > (prog * (nstart/10))) {
          prog <- prog + 1
        }
        setTxtProgressBar(pb, prog * 0.1)
        
        # Starting method
        if(is.null(smartStart)){
          #myseed=seed+run
          #set.seed(myseed)
          set.seed(seed[run])
          randVec= matrix(ceiling(runif(n)*nclus),n,1)
        } else {
          randVec=smartStart
        }
        
        U = tab.disjonctif(randVec)
        # U = data.matrix(fac2disj(randVec))
        #update A
        pseudoinvU = chol2inv(chol(t(U)%*%U))
        P = U%*%pseudoinvU%*%t(U)
        #   R = t(data)%*%((1-alpha)*P-(1-2*alpha)*diag(n))%*%data
        #A = suppressMessages(eigs_sym(R,ndim)$vectors)
        A = eigen(t(data)%*%((1-alpha)*P-(1-2*alpha)*diag(n))%*%data)$vectors[,1:ndim]
        #update Y
        G = data%*%A
        Y = pseudoinvU%*%t(U)%*%G
        f = alpha*ssq(data - G%*%t(A))+(1-alpha)*ssq(data%*%A-U%*%Y)
        f = as.numeric(f) #fixes convergence issue 01 Nov 2016
        fold = f + 2 * conv*f
        iter = 0
        #iterative part
        while (f<fold-conv*f) {
          fold=f
          iter=iter+1
          outK = try(kmeans(G,centers=Y,nstart=100),silent=T)
          
          if(is.list(outK)==FALSE){
            outK = EmptyKmeans(G,centers=Y)
            #  break
          }
          
          v = as.factor(outK$cluster)
          U = diag(nlevels(v))[v,] #dummy cluster membership
          pseudoinvU = chol2inv(chol(t(U)%*%U))
          # update A
          P = U%*%pseudoinvU%*%t(U)
          #R = t(data)%*%((1-alpha)*P-(1-2*alpha)*diag(n))%*%data
          #A = suppressMessages(eigs_sym(R,ndim)$vectors)
          A = eigen(t(data)%*%((1-alpha)*P-(1-2*alpha)*diag(n))%*%data)$vectors
          A = A[,1:ndim]
          G = data %*% A
          #update Y
          Y = pseudoinvU%*%t(U)%*%G
          # criterion
          f = alpha*ssq(data - G%*%t(A))+(1-alpha)*ssq(data%*%A-U%*%Y)
          f = as.numeric(f)
        }
        
        if (f < bestf) {
          bestf = f
          FF = G
          AA = A
          YY = Y
          UU = U
          PP = P
          uu = outK$cluster
        }
      }
      
      cluster = uu##apply(UU, 1, which.max)
      #  size = table(cluster)
      size = table(cluster) 
      aa = sort(size, decreasing = TRUE)
      cluster = mapvalues(cluster, from = as.integer(names(aa)), 
                          to = as.integer(names(table(cluster))))
      centroid = YY#[[mi]]
      centroid = centroid[as.integer(names(aa)), ]
      if (rotation == "varimax") {
        AA = varimax(AA)$loadings[1:m, 1:ndim]
        FF = data%*%AA
        #update center
        centroid = pseudoinverse(t(UU)%*%UU)%*%t(UU)%*%FF 
        centroid = centroid[as.integer(names(aa)),]
      }
      else if (rotation == "promax") {
        AA = promax(AA)$loadings[1:m,1:ndim]
        FF = data%*%AA
        #update center
        centroid = pseudoinverse(t(UU)%*%UU)%*%t(UU)%*%FF
        centroid = centroid[as.integer(names(aa)),]
      }
      
      ##########################
      setTxtProgressBar(pb, 1)
      
      #assign output
      out=list()
      out$obscoord = apply(FF,2, as.numeric) #fixed complex output 16-04-2018
      rownames(out$obscoord) = rownames(data)
      AA = data.matrix(AA)
      out$attcoord = data.matrix(apply(AA,2, as.numeric))#[1:m,1:ndim] 
      rownames(out$attcoord) = colnames(data)
      centroid = data.matrix(centroid)
      out$centroid = apply(centroid, 2, as.numeric) #YY[[mi]]
      names(cluster) = rownames(data)
      out$cluster = cluster #apply(U,1,which.max)
      
      #      Xb = sum(diag((t(A)%*%t(data)%*%U%*%pseudoinverse(crossprod(U))%*%t(U)%*%data%*%A)))
      #      up =  Xb / ((nclus - 1)*ndim + (m-ndim)*ndim)
      #      down = (sum(diag(crossprod(data))) - Xb) / (m*n - nclus*ndim - (m-ndim)*ndim)
      #      out$ch = up / down
      
      out$criterion = bestf
      out$size = as.integer(aa) #round((table(cluster)/sum(table(cluster)))*100,digits=1)
      out$odata = data.frame(odata, stringsAsFactors = TRUE)
      out$scale = scale
      out$center = center
      out$nstart = nstart
      class(out) = "cluspca"
      return(out)
    }
  } else { # dplyr-based implementation
    
    #### A single cluster gives the PCA solution
    if (nclus == 1) {
      nstart = 1
      data = data.frame(data, stringsAsFactors = TRUE)
      n = nrow(data)
      #asymmetric map, biplot
      outp = princomp(data, scale = scale, center = center)
      out=list()
      out$obscoord=outp$scores[,1:ndim] # observations coordinates
      out$attcoord=data.matrix(outp$loadings[,1:ndim]) # attributes coordinates
      rownames(out$obscoord) = rownames(data)
      rownames(out$attcoord) = colnames(data)
      
      out$centroid = 0 #center
      out$cluster = rep(1,n)#cluster
      names(out$cluster) = rownames(data)
      out$criterion = 1 # criterion
      out$size=n #as.integer(aa)  #round((table(cluster)/sum( table(cluster)))*100,digits=1)
      out$odata=data.frame(lapply(data.frame(data, stringsAsFactors = TRUE),factor),stringsAsFactors = TRUE)
      out$scale = scale
      out$center = center
      out$nstart = nstart
      class(out)="cluspca"
      return(out)
    } else {
      #NOTE: FactorialKM needs smartstart k-means or else to perform well
      #FIX: K=2, d=2 does not work for RKM
      if (missing(ndim)) {
        warning('The ndim argument is missing. ndim was set to nclus - 1')
        ndim = nclus - 1
      }
      
      # if (ndim >= nclus) {
      #    stop('The number of clusters should be larger than the number of dimensions.')
      #  }
      
      method <- match.arg(method, c("RKM", "rkm","rKM","FKM", "fkm","fKM"), several.ok = TRUE)[1]
      method <- toupper(method)
      
      if(!is.null(seed)) set.seed(seed)
      seed <- round(2^31 * runif(nstart, -1, 1))
      
      #  If alpha = .5 gives RKM, alpha=1 PCA and alpha =0  FKM.
      if (is.null(alpha) == TRUE)
      {
        if (method == "RKM") {
          alpha = .5
        } else if (method == "FKM") {
          alpha = 0
        }
      }
      odata = data
      data =  scale(data, center = center, scale = scale)
      
      data = data.matrix(data)
      n = dim(data)[1]
      m = dim(data)[2]
      conv=1e-6  # convergence criterion
      bestf = 10^12
      pb <- txtProgressBar(style = 3)
      prog = 1
      
      #    func={}; AA = {}; FF = {}; YY = {}; UU={}
      for (run in c(1:nstart)) {
        if (run > (prog * (nstart/10))) {
          prog <- prog + 1
        }
        setTxtProgressBar(pb, prog * 0.1)
        
        # Starting method
        if(is.null(smartStart)){
          #myseed=seed+run
          #set.seed(myseed)
          set.seed(seed[run])
          randVec= matrix(ceiling(runif(n)*nclus),n,1)
        }else{
          randVec=smartStart
        }
        
        # R = t(data)%*%((1-alpha)*P-(1-2*alpha)*diag(n))%*%data
        # split R = R1 - R2 (12.06.19)
        
        #    U = dummy(randVec)
        #    pseudoinvU = chol2inv(chol(crossprod(U)))
        
        mydata = suppressMessages(as_tibble(cbind(data,group = as.factor(randVec)),.name_repair = "unique"))
        all_groups=tibble(group=mydata$group,trueOrd=1:nrow(mydata))
        # mydata=as_tibble(mydata)
        
        gmeans=mydata%>%
          group_by(group) %>%
          summarise_all(mean)%>%
          full_join(all_groups,gmeans,by="group")%>%
          arrange(trueOrd)%>%
          select(-group,-trueOrd)%>%
          t(.)
        
        R = (1-alpha)*(gmeans)%*%as.matrix(data)
        if (alpha != 0.5) {
          R2 = (1-2*alpha)*crossprod(data)
          R = R - R2
        }
        #A = suppressMessages(eigs(R)$vectors)
        
        #gets ndim + 20% of all dims 
        nd = ndim+round(m*0.2)
        if (nd > ndim) nd = ndim
        if (ncol(R) > 2) 
          A = eigs_sym(R,nd)$vectors
        else
          A =  eigen(R)$vectors
        #    A = eigen(R,symmetric = TRUE)$vectors
        A = A[,1:ndim]
        #update Y
        G = data%*%A
        #  Y = pseudoinvU%*%t(U)%*%G
        all_groups=tibble(group=mydata$group,trueOrd=1:nrow(mydata))
        
        G = suppressMessages(as_tibble(cbind(G,group = as.factor(randVec)),.name_repair = "unique"))
        Y = G%>%
          group_by(group) %>%
          summarise_all(mean) #%>%
        
        
        UY = Y %>%
          full_join(all_groups,Y,by="group")%>%
          arrange(trueOrd)%>%
          select(-group,-trueOrd) #%>%
        
        G = as.matrix(select(G,-group))
        Y = as.matrix(select(Y,-group))
        
        f = alpha*ssq(data - G%*%t(A))+(1-alpha)*ssq(G-UY)
        f = as.numeric(f) #fixes convergence issue 01 Nov 2016
        fold = f + 2 * conv*f
        iter = 0
        #iterative part
        while (f<fold-conv*f) {
          fold=f
          iter=iter+1
          outK = try(kmeans(G,centers=Y,nstart=100),silent=T)
          
          if(is.list(outK)==FALSE){
            outK = EmptyKmeans(G,centers=Y)
            #  break
          }
          
          #  v = as.factor(outK$cluster)
          #  U = diag(nlevels(v))[v,] #dummy cluster membership
          
          #        pseudoinvU = chol2inv(chol(crossprod(U)))
          #  pseudoinvU = chol2inv(chol(t(U)%*%U))
          
          # update A
          mydata = suppressMessages(as_tibble(cbind(data,group = as.factor(outK$cluster)),.name_repair = "unique"))
          all_groups=tibble(group=mydata$group,trueOrd=1:nrow(mydata))
          # mydata=as_tibble(mydata)
          
          gmeans=mydata%>%
            group_by(group) %>%
            summarise_all(mean)%>%
            full_join(all_groups,gmeans,by="group")%>%
            arrange(trueOrd)%>%
            select(-group,-trueOrd)%>%
            t(.)
          R = (1-alpha)*(gmeans)%*%as.matrix(data)
          
          if (alpha != 0.5) {
            R2 = (1-2*alpha)*crossprod(data)
            R = R - R2
          }
          
          #  A = suppressMessages(eigs(R,ncol(data))$vectors)
          if (ncol(R) > 2) 
            A = eigs_sym(R,nd)$vectors
          else
            A =  eigen(R)$vectors
          
          #     A = eigen(R,symmetric = TRUE)$vectors
          A = A[,1:ndim]
          G = data %*% A
          #update Y
          # Y = pseudoinvU%*%t(U)%*%G
          
          #      all_groups=tibble(group=mydata$group,trueOrd=1:nrow(mydata))
          
          G = suppressMessages(as_tibble(cbind(G,group = as.factor(outK$cluster)),.name_repair = "unique"))
          Y = G%>%
            group_by(group) %>%
            summarise_all(mean) #%>%
          
          
          UY = Y %>%
            full_join(all_groups,Y,by="group")%>%
            arrange(trueOrd)%>%
            select(-group,-trueOrd) #%>%
          
          G = as.matrix(select(G,-group))
          Y = as.matrix(select(Y,-group))
          
          # criterion
          f = alpha*ssq(data - G%*%t(A))+(1-alpha)*ssq(G-UY)
          f = as.numeric(f)
        }
        #    func[run] = f
        #    FF[[run]] = G
        #    AA[[run]] = A
        #    YY[[run]] = Y
        #    UU[[run]] = dummy(outK$cluster)
        
        if (f < bestf) {
          bestf = f
          FF = G
          AA = A
          YY = Y
          uu = outK$cluster
        }
        
      }
      
      ##reorder according to cluster size
      UU = tab.disjonctif(uu)
      
      #mi = which.min(func)
      U= UU#UU[[mi]]
      cluster = apply(U,1,which.max)
      
      #csize = round((table(cluster)/sum( table(cluster)))*100,digits=2)
      size = table(cluster)
      aa = sort(size,decreasing = TRUE)
      cluster = mapvalues(cluster, from = as.integer(names(aa)), to = as.integer(names(table(cluster))))
      #reorder centroids
      centroid = YY #YY[[mi]]
      centroid = centroid[as.integer(names(aa)),]
      #######################
      
      ### rotation options ###
      if (rotation == "varimax") { #with Kaiser Normalization
        AA = varimax(AA)$loadings[1:m,1:ndim]
        FF = data%*%AA
        #update center
        centroid = pseudoinverse(t(UU)%*%UU)%*%t(UU)%*%FF 
        centroid = centroid[as.integer(names(aa)),]
        
      } else if (rotation == "promax") {
        AA = promax(AA)$loadings[1:m,1:ndim]
        FF = data%*%AA
        #update center
        centroid = pseudoinverse(t(UU)%*%UU)%*%t(UU)%*%FF
        centroid = centroid[as.integer(names(aa)),]
      }
      
      #  distB = sum(diag(t(AA[[mi]])%*% AA[[mi]]))
      #  distG = sum(diag(t(centroid)%*% centroid))
      #  gamma = ((nclus/m)* distB/distG)^.25
      
      #  AA[[mi]] = (1/gamma)*AA[[mi]]
      #  centroid = gamma*centroid
      #  FF[[mi]] = gamma*FF[[mi]]
      
      ##########################
      setTxtProgressBar(pb, 1)

      #assign output
      out=list()
      out$obscoord = apply(FF,2, as.numeric) #fixed complex output 16-04-2018
      rownames(out$obscoord) = rownames(data)
      AA = data.matrix(AA)
      out$attcoord = data.matrix(apply(AA,2, as.numeric))#[1:m,1:ndim] 
      rownames(out$attcoord) = colnames(data)
      centroid = data.matrix(centroid)
      out$centroid = apply(centroid, 2, as.numeric) #YY[[mi]]
      names(cluster) = rownames(data)
      out$cluster = cluster #apply(U,1,which.max)
      out$criterion = bestf
      out$size = as.integer(aa) #round((table(cluster)/sum(table(cluster)))*100,digits=1)
      out$odata = data.frame(odata,stringsAsFactors = TRUE)
      out$scale = scale
      out$center = center
      out$nstart = nstart
      class(out) = "cluspca"
      return(out)
    }
  }
  
}


ssq = function(a) {
  crossprod(c(as.matrix(a)))
  # t(as.vector(c(as.matrix(a))))%*%as.vector(c(as.matrix(a)))
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

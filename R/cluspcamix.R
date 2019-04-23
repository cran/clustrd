cluspcamix <- function(data, nclus, ndim, method=c("mixedRKM","mixedFKM"), center = TRUE, scale = TRUE, alpha=NULL, rotation="none", nstart=100, smartStart=NULL, seed=NULL, binary = FALSE)
{
  clu=group=trueOrd={}
  . = NULL
  odata = data
  #define X case kapws alliws
  if (binary == FALSE) {
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
    numobs = nrow(data)
    
    #standardize continuous
    # data = scale(data center = center, scale = scale)
    if (center == TRUE) {
      QuantiAct <- t(t(QuantiAct) - as.vector(crossprod(rep(1,numobs)/numobs, as.matrix(QuantiAct))))
    }
    
    if (scale == TRUE) {
      QuantiAct <- t(t(QuantiAct)/sqrt(as.vector(crossprod(rep(1,numobs)/numobs, 
                                                           as.matrix(QuantiAct^2)))))
    }
    lab1a=names(data[, facAct, drop = FALSE])
    lab1b=lapply(data[, facAct, drop = FALSE],function(z) levels(z))
    lab1=rep(lab1a,times=unlist(lapply(lab1b,length)))
    lab2=unlist(lab1b)
    qualilabs=paste(lab1,lab2,sep=".")
    
    QualiAct <-  tab.disjonctif(data[, facAct, drop = FALSE])
    
    
    attlabs = c(colnames(QuantiAct),qualilabs)
    
    
    #standardize categorical
    prop <- colSums(QualiAct * (rep(1,numobs)/numobs))
    # print(as.vector(crossprod(rep(1,numobs)/numobs, as.matrix(QualiAct))))
    
    QualiAct <- t(t(QualiAct) - as.vector(crossprod(rep(1,numobs)/numobs, as.matrix(QualiAct)))) #this is centering MZ
    QualiAct <- t(t(QualiAct)/sqrt(prop))
    
    X <- cbind(QuantiAct, QualiAct)   
  }
  else { #binary = TRUE
    #  numvars <- sapply(data, is.numeric)
    #  anynum <- any(numvars)
    #  catvars <- sapply(data, is.factor)
    #  anyfact <- any(catvars)
    #  if (!anynum) 
    #    cat("\nNo continuous (numeric) variables in data! \n")
    #  if (!anyfact) 
    #    cat("\nNo categorical (factor) variables in data! \n")
    if (is.null(rownames(data))) 
      rownames(data) = 1:nrow(data)
    if (is.null(colnames(data))) 
      colnames(data) = paste("V", 1:ncol(data), sep = "")
    data <- as.data.frame(data)
    data <- droplevels(data)
    numAct <- which(sapply(data, is.numeric))
    facAct <- which(!sapply(data, is.numeric))
    
    QuantiAct <- as.matrix(data[, numAct, drop = FALSE])
    numobs = nrow(data)
    #standardize continuous
    if (center == TRUE) {
      QuantiAct <- t(t(QuantiAct) - as.vector(crossprod(rep(1,numobs)/numobs, as.matrix(QuantiAct))))
    }
    if (scale == TRUE) {
      QuantiAct <- t(t(QuantiAct)/sqrt(as.vector(crossprod(rep(1,numobs)/numobs, 
                                                           as.matrix(QuantiAct^2)))))
    }
    QualiAct <-  data.matrix(data[, facAct, drop = FALSE])#tab.disjonctif(data[, facAct, drop = FALSE])
    
    #  lab1a=names(data[, facAct, drop = FALSE])
    #  lab1b=lapply(data[, facAct, drop = FALSE],function(z) levels(z))
    #  lab1=rep(lab1a,times=unlist(lapply(lab1b,length)))
    #  lab2=unlist(lab1b)
    #  qualilabs=paste(lab1,lab2,sep=".")
    
    attlabs = c(colnames(QuantiAct),colnames(QualiAct))
    #standardize categorical
    prop <- colSums(QualiAct * (rep(1,numobs)/numobs))
    
    QualiAct <- t(t(QualiAct) - as.vector(crossprod(rep(1,numobs)/numobs, as.matrix(QualiAct)))  ) #this is centering MZ
    QualiAct <- t(t(QualiAct)/sqrt(prop))
    
    X <- cbind(QuantiAct, QualiAct)   
    #   str(head(data.frame(X)))
  }
  
  if(!is.null(seed)) set.seed(seed)
  seed <- round(2^31 * runif(nstart, -1, 1))
  pb <- txtProgressBar(style = 3)
  prog = 1
  
  data = X
  
  if (nrow(data) < 700) { #threshold for largish data sets
    
    #######
    if (nclus == 1) {
      nstart = 1
      data = data.frame(data)
      n = nrow(data)
      outp = princomp(data)
      out=list()
      out$obscoord = outp$scores[,1:ndim] # observations coordinates
      out$attcoord = data.matrix(outp$loadings[,1:ndim]) # attributes coordinates
      rownames(out$obscoord) = rownames(data)
      rownames(out$attcoord) = attlabs#colnames(data)
      
      out$centroid = 0 #center
      out$cluster = rep(1,n)#cluster
      names(out$cluster) = rownames(data)
      out$criterion = 1 # criterion
      out$size=n #as.integer(aa)  #round((table(cluster)/sum( table(cluster)))*100,digits=1)
      out$odata = odata #data.frame(data) #data.frame(lapply(data.frame(data),factor))
      #   out$scale = FALSE
      #    out$center = FALSE
      out$nstart = nstart
      class(out) = "cluspcamix"
      return(out)
    } else {
      
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
      
      method <- toupper(method)
      
      method <- match.arg(method, c("MIXEDRKM", "MIXEDFKM"), several.ok = T)[1]
      
      
      if (is.null(alpha) == TRUE) {
        if (method == "MIXEDRKM") {
          alpha = 0.5
        }
        else if (method == "MIXEDFKM") 
          alpha = 0
      }
      
      # data = scale(data center = center, scale = scale)
      data = data.matrix(data)
      n = nrow(data) #dim(data)[1]
      m = ncol(data) #dim(data)[2]
      conv = 1e-06
      bestf = 10^12
      pb <- txtProgressBar(style = 3)
      prog = 1
      for (run in c(1:nstart)) {
        if (run > (prog * (run/9))) {
          prog <- prog + 1
        }
        setTxtProgressBar(pb, prog * 0.1)
        
        
        if (is.null(smartStart)) {
          #  myseed = seed + run
          #  set.seed(myseed)
          set.seed(seed[run])
          randVec = matrix(ceiling(runif(n) * nclus), n, 
                           1)
        }
        else {
          randVec = smartStart
        }
        #U = dummy(randVec)
        U = data.matrix(fac2disj(randVec))
        
        P = U %*% pseudoinverse(t(U) %*% U) %*% t(U)
        A = eigen(t(data) %*% ((1 - alpha) * P - (1 - 2 * 
                                                    alpha) * diag(n)) %*% data)$vectors
        A = A[, 1:ndim]
        G = data %*% A
        Y = pseudoinverse(t(U) %*% U) %*% t(U) %*% G
        f = alpha * ssq(data - G %*% t(A)) + (1 - alpha) * 
          ssq(data %*% A - U %*% Y)
        f = as.numeric(f)
        fold = f + 2 * conv * f
        iter = 0
        while (f < fold - conv * f) {
          fold = f
          iter = iter + 1
          outK = try(kmeans(G, centers = Y, nstart = 100), 
                     silent = T)
          if (is.list(outK) == FALSE) {
            outK = EmptyKmeans(G, centers = Y)
          }
          v = as.factor(outK$cluster)
          U = diag(nlevels(v))[v, ]
          P = U %*% pseudoinverse(t(U) %*% U) %*% t(U)
          A = eigen(t(data) %*% ((1 - alpha) * P - (1 - 
                                                      2 * alpha) * diag(n)) %*% data)$vectors
          A = A[, 1:ndim]
          G = data %*% A
          Y = pseudoinverse(t(U) %*% U) %*% t(U) %*% G
          f = alpha * ssq(data - G %*% t(A)) + (1 - alpha) * 
            ssq(data %*% A - U %*% Y)
          f = as.numeric(f)
        }
        if (f < bestf) {
          bestf = f
          FF = G
          AA = A
          YY = Y
          uu = outK$cluster
        }
        
      }
      #  mi = which.min(func)
      ##reorder according to cluster size
      UU = dummy(uu)
      
      U = UU#[[mi]]
      cluster = apply(U, 1, which.max)
      size = table(cluster)
      aa = sort(size, decreasing = TRUE)
      cluster = mapvalues(cluster, from = as.integer(names(aa)), 
                          to = as.integer(names(table(cluster))))
      centroid = YY#[[mi]]
      centroid = centroid[as.integer(names(aa)), ]
      if (rotation == "varimax") {
        AA = varimax(AA)$loadings[1:m, 1:ndim]
        FF = data %*% AA
        centroid = pseudoinverse(t(U) %*% U) %*% t(U) %*% 
          FF
        centroid = centroid[as.integer(names(aa)), ]
      }
      else if (rotation == "promax") {
        AA = promax(AA)$loadings[1:m, 1:ndim]
        FF = data %*% AA
        centroid = pseudoinverse(t(U) %*% U) %*% t(U) %*% 
          FF
        centroid = centroid[as.integer(names(aa)), ]
      }
      setTxtProgressBar(pb, 1)
      
      out = list()
      #  mi = which.min(func)
      out$obscoord = FF
      rownames(out$obscoord) = rownames(data)
      out$attcoord = data.matrix(AA)
      rownames(out$attcoord) = attlabs#colnames(data)
      out$centroid = centroid
      names(cluster) = rownames(data)
      out$cluster = cluster
      out$criterion = bestf #func[mi]
      out$size = as.integer(aa)
      out$odata = odata #data.frame(data) #check what's needed for tunecluspcamix
      # colnames(out$odata) = #rownames(out$attcoord)
      #print(colnames(out$odata))
      out$scale = scale
      out$center = center
      out$nstart = nstart
      class(out) = "cluspcamix"
      return(out)
    }
  } else { # dplyr-based implementation
    
    #######
    if (nclus == 1) {
      nstart = 1
      data = data.frame(data)
      n = nrow(data)
      outp = princomp(data)
      out=list()
      out$obscoord = outp$scores[,1:ndim] # observations coordinates
      out$attcoord = data.matrix(outp$loadings[,1:ndim]) # attributes coordinates
      rownames(out$obscoord) = rownames(data)
      rownames(out$attcoord) = attlabs#colnames(data)
      
      out$centroid = 0 #center
      out$cluster = rep(1,n)#cluster
      names(out$cluster) = rownames(data)
      out$criterion = 1 # criterion
      out$size=n #as.integer(aa)  #round((table(cluster)/sum( table(cluster)))*100,digits=1)
      out$odata = odata #data.frame(data) #data.frame(lapply(data.frame(data),factor))
      #   out$scale = FALSE
      #    out$center = FALSE
      out$nstart = nstart
      class(out) = "cluspcamix"
      return(out)
    } else {
      
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
      
      method <- toupper(method)
      
      method <- match.arg(method, c("MIXEDRKM", "MIXEDFKM"), several.ok = T)[1]
      
      
      if (is.null(alpha) == TRUE) {
        if (method == "MIXEDRKM") {
          alpha = 0.5
        }
        else if (method == "MIXEDFKM") 
          alpha = 0
      }
      
      # data = scale(data center = center, scale = scale)
      data = data.matrix(data)
      n = nrow(data) #dim(data)[1]
      m = ncol(data) #dim(data)[2]
      conv = 1e-06
      bestf = 10^12
      pb <- txtProgressBar(style = 3)
      prog = 1
      for (run in c(1:nstart)) {
        if (run > (prog * (run/9))) {
          prog <- prog + 1
        }
        setTxtProgressBar(pb, prog * 0.1)
        
        
        if (is.null(smartStart)) {
          #  myseed = seed + run
          #  set.seed(myseed)
          set.seed(seed[run])
          randVec = matrix(ceiling(runif(n) * nclus), n, 
                           1)
        }
        else {
          randVec = smartStart
        }
       
        mydata = as_tibble(cbind(data,group = as.factor(randVec)))
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
        #A = suppressWarnings(eigs(R)$vectors)
        
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
        
        G = as_tibble(cbind(G,group = as.factor(randVec)))
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
        
        fold = f + 2 * conv * f
        iter = 0
        while (f < fold - conv * f) {
          fold = f
          iter = iter + 1
          outK = try(kmeans(G, centers = Y, nstart = 100), 
                     silent = T)
          if (is.list(outK) == FALSE) {
            outK = EmptyKmeans(G, centers = Y)
          }
          v = as.factor(outK$cluster)
          U = diag(nlevels(v))[v, ]
          P = U %*% pseudoinverse(t(U) %*% U) %*% t(U)
          A = eigen(t(data) %*% ((1 - alpha) * P - (1 - 
                                                      2 * alpha) * diag(n)) %*% data)$vectors
          A = A[, 1:ndim]
          G = data %*% A
          Y = pseudoinverse(t(U) %*% U) %*% t(U) %*% G
          f = alpha * ssq(data - G %*% t(A)) + (1 - alpha) * 
            ssq(data %*% A - U %*% Y)
          f = as.numeric(f)
        }
        if (f < bestf) {
          bestf = f
          FF = G
          AA = A
          YY = Y
          uu = outK$cluster
        }
        
      }
      #  mi = which.min(func)
      ##reorder according to cluster size
      UU = dummy(uu)
      
      U = UU#[[mi]]
      cluster = apply(U, 1, which.max)
      size = table(cluster)
      aa = sort(size, decreasing = TRUE)
      cluster = mapvalues(cluster, from = as.integer(names(aa)), 
                          to = as.integer(names(table(cluster))))
      centroid = YY#[[mi]]
      centroid = centroid[as.integer(names(aa)), ]
      if (rotation == "varimax") {
        AA = varimax(AA)$loadings[1:m, 1:ndim]
        FF = data %*% AA
        centroid = pseudoinverse(t(U) %*% U) %*% t(U) %*% 
          FF
        centroid = centroid[as.integer(names(aa)), ]
      }
      else if (rotation == "promax") {
        AA = promax(AA)$loadings[1:m, 1:ndim]
        FF = data %*% AA
        centroid = pseudoinverse(t(U) %*% U) %*% t(U) %*% 
          FF
        centroid = centroid[as.integer(names(aa)), ]
      }
      setTxtProgressBar(pb, 1)
      
      out = list()
      #  mi = which.min(func)
      out$obscoord = FF
      rownames(out$obscoord) = rownames(data)
      out$attcoord = data.matrix(AA)
      rownames(out$attcoord) = attlabs#colnames(data)
      out$centroid = centroid
      names(cluster) = rownames(data)
      out$cluster = cluster
      out$criterion = bestf #func[mi]
      out$size = as.integer(aa)
      out$odata = odata #data.frame(data) #check what's needed for tunecluspcamix
      # colnames(out$odata) = #rownames(out$attcoord)
      #print(colnames(out$odata))
      out$scale = scale
      out$center = center
      out$nstart = nstart
      class(out) = "cluspcamix"
      return(out)
    }
  }
  
  
}

moy.ptab <- function(V, poids) {
  as.vector(crossprod(poids/sum(poids), as.matrix(V)))
}

ec.tab <- function(V, poids) {
  ecart.type <- sqrt(as.vector(crossprod(poids/sum(poids), 
                                         as.matrix(V^2))))
  ecart.type[ecart.type <= 1e-16] <- 1
  return(ecart.type)
}

fac2disj<- function(fac, drop = FALSE) {
  ## Returns the disjunctive table corrseponding to a factor
  n <- length(fac)
  fac <- as.factor(fac)
  if(drop)
    fac <- factor(fac)
  x <- matrix(0, n, nlevels(fac))
  x[(1:n) + n * (unclass(fac) - 1)] <- 1
  dimnames(x) <- list(names(fac), as.character(levels(fac)))
  return(data.frame(x, check.names = FALSE))
}

tab.disjonctif.NA <-
  function (tab,rename.level=FALSE) {
    tab <- as.data.frame(tab)
    modalite.disjonctif <- function(i) {
      moda <- tab[, i]
      nom <- names(tab)[i]
      n <- length(moda)
      moda <- as.factor(moda)
      x <- matrix(0, n, length(levels(moda)))
      ind<-(1:n) + n * (unclass(moda) - 1)
      indNA<-which(is.na(ind))
      
      x[(1:n) + n * (unclass(moda) - 1)] <- 1
      x[indNA,]<-NA 
      if (rename.level==TRUE){
        dimnames(x) <- list(row.names(tab), paste(nom, levels(moda), 
                                                  sep = "="))
      }    else{
        dimnames(x) <- list(row.names(tab), levels(moda))
      }
      
      
      return(x)
    }
    if (ncol(tab) == 1) 
      res <- modalite.disjonctif(1)
    else {
      res <- lapply(1:ncol(tab), modalite.disjonctif)
      res <- as.matrix(data.frame(res, check.names = FALSE))
    }
    return(res)
  }

ssq = function(a) {
  t(as.vector(c(as.matrix(a))))%*%as.vector(c(as.matrix(a)))
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

tab.disjonctif <- function(tab) {
  tab <- as.data.frame(tab)
  modalite.disjonctif <- function(i) {
    moda <- as.factor(tab[, i])
    n <- length(moda)
    x <- matrix(0L, n, nlevels(moda))
    x[(1:n) + n * (unclass(moda) - 1L)] <- 1L
    return(x)
  }
  if (ncol(tab) == 1) {
    res <- modalite.disjonctif(1)
    dimnames(res) <- list(attributes(tab)$row.names, levels(tab[, 
                                                                1]))
  }
  else {
    variable <- rep(attributes(tab)$names, lapply(tab, nlevels))
    listModa <- unlist(lapply(tab, levels))
    wlistModa <- which((listModa) %in% c("y", "n", "Y", "N"))
    if (!is.null(wlistModa)) 
      listModa[wlistModa] <- paste(variable[wlistModa], 
                                   listModa[wlistModa], sep = ".")
    numlistModa <- which(unlist(lapply(listModa, is.numeric)))
    if (!is.null(numlistModa)) 
      listModa[numlistModa] <- paste(variable[numlistModa], 
                                     listModa[numlistModa], sep = ".")
    res <- lapply(1:ncol(tab), modalite.disjonctif)
    res <- as.matrix(data.frame(res, check.names = FALSE))
    dimnames(res) <- list(attributes(tab)$row.names, listModa)
  }
  return(res)
}

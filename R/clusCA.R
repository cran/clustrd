clusCA <- function(data,nclus,ndim,nstart=100,smartStart=NULL,gamma = FALSE, seed=NULL, inboot= FALSE){
  K = nclus
  k = ndim
  nrs = nstart
  q = ncol(data)
  maxiter = 100
  maxinert=-1
  data = data.frame(data)
  if (inboot == FALSE) {
    data=as.data.frame(lapply(data,as.factor))
    Z = tab.disjonctif(data)#dummy.data.frame(data, dummy.classes = "ALL") # The original super indicator
    lab1a=names(data)
    lab1b=lapply(data,function(z) levels(as.factor(z)))
    lab1=rep(lab1a,times=unlist(lapply(lab1b,length)))
    lab2=unlist(lab1b)
    colnames(Z) = paste(lab1,lab2,sep=".")
  }
  else
    Z = data
  n=nrow(Z)
  Q=ncol(Z)
  
  Dzh=diag(as.vector(colSums(Z)^(.5))) 
  Dzhi=pseudoinverse(Dzh)
  #Dzhi= chol2inv(chol(Dzh))
  MZ=scale(Z,scale=FALSE)
  MZD = MZ %*% Dzhi
  if(!is.null(seed)) set.seed(seed)
  seed <- round(2^31 * runif(nstart, -1, 1))
  
  pb <- txtProgressBar(style = 3)
  prog = 1
  # Do nrs random starts
  fvec=c()
  for (rs in 1:nrs){
    if (rs > (prog * (nrs/10))) {
      prog <- prog + 1
    }
    setTxtProgressBar(pb, prog * 0.1)
    
    if(is.null(smartStart)){
      # myseed=seed+rs
      #  set.seed(myseed)
      set.seed(seed[rs])
      randVec= matrix(ceiling(runif(n)*nclus),n,1)
    }else{
      randVec=smartStart
    }
    
   Zki=tab.disjonctif(randVec)

    Dk=t(Zki) %*% Zki
    Dks=Dk^(.5)
    #Dksi=pseudoinverse(Dks)
    Dksi= chol2inv(chol(Dks))
    DZkZD= sqrt(n/q)*Dksi%*% t(Zki) %*% MZD     # equation 6 on the paper
    svdDZkZD=svd(DZkZD)
    #this is for ndim = 1 to work
    if (k != 1) {
      Lk=diag(svdDZkZD$d[1:k])
    } else {
      Lk = data.matrix(svdDZkZD$d[1])
    }
    G = svdDZkZD$u

    Gi = G[,1:k]

    Gi=Dksi%*%Gi%*%Lk # CA row coordinates (section 2 in the paper right below formula (1))
  #  Gs = Gi%*%pseudoinverse(Lk)
  #  ARC = Dk%*%(Gs*Gs)# cluster means absolute contributions
    
    Bstar=svdDZkZD$v
    B=sqrt(n*q)*Dzhi %*% Bstar # as in eq. (8)
    
    Bns=B[,1:k]          
    Bi=Bstar[,1:k]           # attribute quantifications. Orthonormal
    
    #              
    #  inertia = sum(t(Lk)%*% Lk) # explained inertia in k dimensions
    
    Yi=sqrt((n/q)) * MZD%*%Bi  # The coordinates for the subjects. as in eq. (10).
 #   Ys = Zki%*%chol2inv(chol(Dk))%*%t(Zki)%*%Yi%*%pseudoinverse(Lk)
#    ARCy = Ys*Ys
    
    #######
#    ARCb = Bi*Bi
#    print(ARCb)
#    print(apply(ARCb,2,sum))
    #######
    
    GDGbef=sum(diag(t(Gi) %*% Dk %*% Gi))   # Objective value before K means This is equivalent to formula (6), but then in the paper
    objbef=GDGbef
    ######## END first fixed C step: Given random C, B and G are optimal
    improv=10
    iter=0   
    
  #  objective=sum(diag((t(Yi)%*%Zki%*%chol2inv(chol(t(Zki)%*%Zki))%*%t(Zki)%*%Yi)))
    while ((improv > 0.0001) && (iter<=maxiter)){
      iter = iter + 1
      outK = try(kmeans(Yi,centers=Gi,nstart=100),silent=T)
      #empty clusters
      if(is.list(outK) == F){
        outK=EmptyKmeans(Yi,centers=Gi) 
        #  break 
      }
      Zki=outK$cluster
      Gi=outK$centers
      Zki = tab.disjonctif(Zki)
      
      Dk = t(Zki) %*% Zki        
      Dks = Dk^(.5)                  # New Dch weights
      #Dksi = pseudoinverse(Dks)
      Dksi = chol2inv(chol(Dks))
      #END Fixed B step: Given B, C is optimal.
      #Now: Fix C and recalculate B
      DkZkZD=sqrt(n/q)*Dksi%*% t(Zki)%*%MZD # Make the new matrix for the SVD.
      
      outDkZkZD=svd(DkZkZD)       # New SVD, CA analysis
      G=outDkZkZD$u
      Gi=G[,1:k]
      # this is for ndim = 1 to work
      if (k != 1) {
        Lk=diag(outDkZkZD$d[1:k])
      } else {
        Lk = data.matrix(outDkZkZD$d[1])
      }
      
      Gi = Dksi %*% Gi %*% Lk
      
#      Gs = Gi%*%pseudoinverse(Lk)
#      ARC = Dk%*%(Gs*Gs)# cluster means absolute contributions
      
      Bstar=outDkZkZD$v
      Bi=Bstar[,1:k]
      B=sqrt(n*q)*Dzhi%*%Bi
      Bns=B[,1:k]
 
      Yi= sqrt((n/q)) * MZD %*% Bi # Subject coordinates
#      Ys = Zki%*%chol2inv(chol(Dk))%*%t(Zki)%*%Yi%*%pseudoinverse(Lk)
#      ARCy = Ys*Ys
      
      improv=sum(diag(t(Gi)%*% t(Zki)%*%Zki%*%Gi))-objbef
      objbef=sum(diag(t(Gi) %*% t(Zki) %*% Zki %*% Gi))
      #    print(improv)
      varsi=sum(diag(t(Gi)%*% Dk %*% Gi)) #no need to calc inside the loop
      
    }
    #FIX: 16/09/2016, calculated outside the loop
    # varsi=sum(diag(t(Gi)%*% Dk %*% Gi))
    
    # fvec = c(fvec, objbef)
    #  print(fvec)
    #  print(rs)
    if (varsi > maxinert){ #gamma
      if (gamma == TRUE) { 
        distB = sum(diag(t(Bns)%*%  Bns))
        distG = sum(diag(t(Gi)%*% Gi))
        g = ((K/Q)* distB/distG)^.25
  #      ARC = ARC
        
  #      ARCy = ARCy
  #      G2 = Dk%*%(G[,1:nclus-1]*G[,1:nclus-1])
  #      totg = apply(G2,1,sum)
  #      RPCG = pseudoinverse(diag(totg))%*%G2
        
  #      Y2 = Yi*Yi
  #      toty = apply(Y2,1,sum)
  #      RRCY = chol2inv(chol(diag(toty)))%*%Y2
        Bsol = (1/g)*Bns
        Gsol = g*Gi
        Ysol = g*Yi
      } else {
   #     ARC = ARC
  #      ARCy = ARCy
        
   #     G2 = Dk%*%(G[,1:nclus-1]*G[,1:nclus-1])
  #      totg = apply(G2,1,sum)
    #    RPCG = pseudoinverse(diag(totg))%*%G2
        
    #    Y2 = Yi*Yi
    #    toty = apply(Y2,1,sum)
    #    RRCY= chol2inv(chol(diag(toty)))%*%Y2
        
        Bsol = Bns
        Gsol = Gi
        Ysol = Yi
      }
      Csol = Zki
      maxinert = varsi
    }
  }
  
  wone = which(Csol==1,arr.ind=T)
  cluster = matrix(0,n,1)
  cluster[wone[,1]] = wone[,2]
  
  ##reorder cluster membership according to cluster size
  #csize = round((table(cluster)/sum( table(cluster)))*100,digits=2)
  #aa = sort(csize,decreasing = TRUE)
  size = table(cluster)
  aa = sort(size,decreasing = TRUE)
  cluster = mapvalues(cluster, from = as.integer(names(aa)), to = as.integer(names(table(cluster))))
  #reorder centroids
  Gsol = Gsol[as.integer(names(aa)),]
  setTxtProgressBar(pb, 1)
  out=list()
  out$obscoord=data.frame(Ysol) # observations coordinates
  out$attcoord=data.frame(Bsol) # attributes coordinates
  out$centroid=data.frame(Gsol) # centroids
 # out$ARC=data.frame(ARC)
#  out$ARCy=data.frame(ARCy)
#  out$RPCG=data.frame(RPCG)
#  out$RRCY=data.frame(RRCY)
  cluster = as.integer(cluster)
  names(cluster) = rownames(data) 
  rownames(out$obscoord) = rownames(data)
  rownames(out$attcoord) = colnames(Z) 
  
  out$cluster=cluster   # cluster membership
  out$criterion=maxinert # criterion
  #  out$iters=iters # number of iterations
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

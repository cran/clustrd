clusCA <- function(data,nclus,ndim,nstart=100,smartStart=NULL,gamma = FALSE, seed=NULL, inboot= FALSE){
  K = nclus
  k = ndim
  nrs = nstart
  q = ncol(data)
  maxiter = 100
  maxinert=-1
  data = data.frame(data, stringsAsFactors = TRUE)
  odata = data
  if (inboot == FALSE) {
    data=as.data.frame(lapply(data,as.factor))
    Z = tab.disjonctif(data)#dummy.data.frame(data, dummy.classes = "ALL") # The original super indicator
    lab1a=names(data)
    lab1b=lapply(data,function(z) levels(as.factor(z)))
    lab1=abbreviate(rep(lab1a,times=unlist(lapply(lab1b,length))),3)
    lab2=unlist(lab1b)
    colnames(Z) = paste(lab1,lab2,sep=".")
  }
  else
    Z = data
  

  n=nrow(Z)
  Q=ncol(Z)
  
  Dz=diag(apply(Z,2,sum)) # Dz
  Dzhi=pseudoinverse(Dz^(0.5)) # %Dz^-.5
  
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
    Dksi= pseudoinverse(Dks)#chol2inv(chol(Dks))
    #sqrt(n/q)
    DZkZD= Dksi%*% t(Zki) %*% MZD     # equation 6 on the paper
    svdDZkZD=svd(DZkZD)
    
    #this is for ndim = 1 to work
    if (k != 1) {
      Lk=diag(svdDZkZD$d[1:k])
    } else {
      Lk = data.matrix(svdDZkZD$d[1])
    }
  
    G = data.matrix(svdDZkZD$u[,1:k])
    Gi = G[,1:k]
    Gi=Dksi%*%Gi%*%Lk # CA row coordinates (section 2 in the paper right below formula (1))
    #inertia=sum(diag(Lk%*%Lk)) # explained inertia in k dimensions
    Bstar = data.matrix(svdDZkZD$v[,1:k])
    #sqrt(n*q)*
    B=Dzhi %*% Bstar # as in eq. (8)
    Bns = B[,1:k]
    Bi=Bstar[,1:k]           # attribute quantifications. Orthonormal
    #sqrt((n/q))*
    Yi= MZD%*%Bi  # The coordinates for the subjects. as in eq. (10).
    
    GDGbef=sum(diag(t(Gi) %*% Dk %*% Gi))   # Objective value before K means This is equivalent to formula (6), but then in the paper
    objbef=GDGbef
    ######## END first fixed C step: Given random C, B and G are optimal
    improv=10
    iter=0   
    
    #  objective=sum(diag((t(Yi)%*%Zki%*%chol2inv(chol(t(Zki)%*%Zki))%*%t(Zki)%*%Yi)))
    while ((improv > 0.0001) && (iter<maxiter)){
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
      Dksi = pseudoinverse(Dks)
      #Dksi = chol2inv(chol(Dks))
      #END Fixed B step: Given B, C is optimal.
      #Now: Fix C and recalculate B
      #
      #sqrt(n/q)
      DkZkZD=Dksi%*% t(Zki)%*%MZD # Make the new matrix for the SVD.
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
      Gs = Gi%*%pseudoinverse(Lk)
  #    abscluscontr = Dk%*%(Gs*Gs)# cluster means absolute contributions
      
      Bstar=outDkZkZD$v
      
      Bi=Bstar[,1:k]
      #sqrt(n*q)
      B=Dzhi%*%Bstar
      Bns=B[,1:k]
  #    absattcontr = Dz%*%Bns*Bns
      
  #    Bp = B[,1:(nclus-1)]%*%diag(outDkZkZD$d)[1:nclus-1,1:nclus-1]
  #    B2= Dzhi%*%(Bp[,1:nclus-1]*Bp[,1:nclus-1]) 
  #    totb= apply(B2,1,sum)
  #    relattcontr = pseudoinverse(diag(totb))%*%B2
      #sqrt((n/q)) *
      Yi=  MZD %*% Bi # Subject coordinates
   #   Ys = Zki%*%chol2inv(chol(Dk))%*%t(Zki)%*%Yi%*%pseudoinverse(Lk)
  #    absobscontr = Ys*Ys
      
      improv=sum(diag(t(Gi)%*% t(Zki)%*%Zki%*%Gi))-objbef
      objbef=sum(diag(t(Gi) %*% t(Zki) %*% Zki %*% Gi))
      varsi=sum(diag(t(Gi)%*% Dk %*% Gi)) #no need to calc inside the loop
      
    }
  
    if (varsi > maxinert){ #gamma
      if (gamma == TRUE) { 
        distB = sum(diag(t(Bns)%*%  Bns))
        distG = sum(diag(t(Gi)%*% Gi))
        g = ((K/Q)* distB/distG)^.25
   #     abscluscontr = abscluscontr
   #    absattcontr = absattcontr
   #     absobscontr = absobscontr
        
    #    G2 = Dk%*%(G[,1:nclus-1]*G[,1:nclus-1])
    #    totg = apply(G2,1,sum)
    #    relcluscontr = pseudoinverse(diag(totg))%*%G2
        
    #    relattcontr = relattcontr
        
    #    Y2 = Yi*Yi
    #    toty = apply(Y2,1,sum)
    #    relobscontr = chol2inv(chol(diag(toty)))%*%Y2
        
        Bsol = (1/g)*Bns
        Gsol = g*Gi
        Ysol = g*Yi
      } else {
    #    abscluscontr = abscluscontr
    #    absobscontr = absobscontr
    #    absattcontr = absattcontr
        
    #    G2 = Dk%*%(G[,1:nclus-1]*G[,1:nclus-1])
    #    totg = apply(G2,1,sum)
    #    relcluscontr = pseudoinverse(diag(totg))%*%G2
        
    #    relattcontr = relattcontr
        
    #    Y2 = Yi*Yi
    #    toty = apply(Y2,1,sum)
    #    relobscontr= chol2inv(chol(diag(toty)))%*%Y2
        
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
  out$obscoord=data.frame(Ysol,stringsAsFactors = TRUE) # observations coordinates
  out$attcoord=data.frame(Bsol,stringsAsFactors = TRUE) # attributes coordinates
  out$centroid=data.frame(Gsol,stringsAsFactors = TRUE) # centroids
#  out$abscluscontr=data.frame(abscluscontr)
#  out$absobscontr=data.frame(absobscontr)
#  out$absattcontr=data.frame(absattcontr)

#  out$relcluscontr=data.frame(relcluscontr)
#  out$relattcontr=data.frame(relattcontr)
#  out$relobscontr=data.frame(relobscontr)
  cluster = as.integer(cluster)
  names(cluster) = rownames(data) 
#  rownames(out$obscoord) = rownames(data)
#  rownames(out$absobscontr) = rownames(data)
#  rownames(out$relobscontr) = rownames(data)
  
#  rownames(out$attcoord) = colnames(Z) 
#  rownames(out$absattcontr) = colnames(Z) 
#  rownames(out$relattcontr) = colnames(Z) 

  out$cluster=cluster   # cluster membership
  out$criterion=maxinert # criterion
  #  out$iters=iters # number of iterations
  out$size=as.integer(aa) #round((table(cluster)/sum( table(cluster)))*100,digits=1)
  out$odata=data.frame(odata, stringsAsFactors = TRUE)#data.frame(lapply(data.frame(data,stringsAsFactors = TRUE),factor),stringsAsFactors = TRUE)
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

EmptyKmeans<-function(data,centers){
  out = list()
  n = nrow(data)
  if(dim(centers)[1]==1){
    K=centers
    print(K)
    start=ceiling(runif(K)*n)
    centers=data[start,]
  }else{
    K=nrow(centers)
    
  }
  gvec=matrix(1,nrow=n,1)
  it=0
  maxiter=1000
  oldf=10000000
  gmat=c()
  while(it<=maxiter){
    it=it+1
    
    for(i in 1:n){
      assignVec=c()
      for(j in 1:K){
        a=data[i,]
        b=centers[j,]
        assignVec[j]=dist(rbind(a,b))
      }
      gvec[i]=which.min(assignVec)
      
    }
    
    C=dummy(gvec)
    
    G=chol2inv(chol(t(C)%*% C))%*% t(C) %*% data

    CG=C%*%G
    if(length( table(gvec))!=K){
      #print("there is an empty cluster")
      who_empty=setdiff(1:K,unique(gvec))
      #print(who_empty)
      allDist=sqrt(apply((data-CG)^2,1,sum))
      for (emp in 1:length(who_empty)){
        far=which.min(allDist)
        while(length(which(gvec==gvec[far]))==1){
          allDist=allDist[-far]
          far=which.min(allDist)
        }
        gvec[far]=who_empty[emp]
        #print(unique(gvec))
        allDist=allDist[-far]
      }
      
      C=dummy(gvec)
      G=chol2inv(chol(t(C)%*% C))%*% t(C) %*% data
      CG=C%*%G
    }
    gmat=cbind(gmat,gvec)
    f = sum(diag(t(data-CG) %*% (data-CG)))
    f = as.numeric(f) #fixes complex value 16-04-2018
    #print(f)
    if((oldf-f) <= 0){
      out$cluster=gvec
      out$centers=G
      out$f=f
      it = maxiter + 1
    }
    oldf = f
    centers = G
  }
  
  out
}
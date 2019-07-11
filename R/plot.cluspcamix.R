plot.cluspcamix<-function(x, dims = c(1,2), cludesc = FALSE, topstdres = 20, attlabs = NULL, subplot = FALSE, what = c(TRUE,TRUE), ...){
  
  d1 = NULL
  d2 = NULL
  gr = NULL
  olab = NULL
  act = NULL
  attnam = NULL
  slp = NULL
  out=list()
  if (dim(data.frame(x$attcoord))[2] == 1) {
    stop('There is only one dimension. A 2D scatterplot cannot be produced.')
  } 
  # out=list()
  dim1=dims[1]
  dim2=dims[2]
  K=max(x$cluster)
  
  if (missing(attlabs)) {
    attlabs=rownames(x$attcoord)#colnames(x$odata) #row.names(x$attcoord)
  }
  
  #do not show obs labels if more than 30
  if (dim(x$obscoord)[1] < 30) {
    obslabs = row.names(x$odata)
  } else
  {
    obslabs = paste("")
  }
  xallmax=max(max(x$attcoord[,dim1]),max(x$obscoord[,dim1]))
  xallmin=min(min(x$attcoord[,dim1]),min(x$obscoord[,dim1]))
  yallmax=max(max(x$attcoord[,dim2]),max(x$obscoord[,dim2]))
  yallmin=min(min(x$attcoord[,dim2]),min(x$obscoord[,dim2]))
  
  #pdf(file=paste("K",deparse(K),"Mapunits.pdf",sep=""),height=9 , width=9)
  
  xallmax=max(max(x$attcoord[,dim1]),max(x$obscoord[,dim1]))
  xallmin=min(min(x$attcoord[,dim1]),min(x$obscoord[,dim1]))
  yallmax=max(max(x$attcoord[,dim2]),max(x$obscoord[,dim2]))
  yallmin=min(min(x$attcoord[,dim2]),min(x$obscoord[,dim2]))
  x_all_range=xallmax-xallmin
  y_all_range=yallmax-yallmin
  all_range=max(x_all_range,y_all_range)
  xallmax=xallmin+all_range
  yallmax=yallmin+all_range
  # extend borders to show all points
  # xallmax = xallmax + 0.05*all_range
  #  xallmin = xallmin - 0.05*all_range
  #  yallmax = yallmax + 0.05*all_range
  #  yallmin = yallmin - 0.05*all_range
  #
  xattmax=max(x$attcoord[,dim1])
  xattmin=min(x$attcoord[,dim1])
  yattmax=max(x$attcoord[,dim2])
  yattmin=min(x$attcoord[,dim2])
  x_att_range=xattmax-xattmin
  y_att_range=yattmax-yattmin
  att_range=max(x_att_range,y_att_range)
  xattmax=xattmin+att_range
  yattmax=yattmin+att_range
  ######################################################
  att_df=data.frame(d1=x$attcoord[,dim1],d2=x$attcoord[,dim2],attnam=attlabs)
  glab=paste(rep("C",K),1:K,sep="")  
  if (length(x$size) != 1)
  {
    group_df= data.frame(d1=x$centroid[,dim1],d2=x$centroid[,dim2],glab=glab,gr=levels(factor(x$cluster)))
  }
  obs_df=data.frame(d1=x$obscoord[,dim1],d2=x$obscoord[,dim2],gr=factor(x$cluster),olab=obslabs)
  xact=union(which(att_df$d1> .5),which(att_df$d1< -.5))
  yact=union(which(att_df$d2>.5), which(att_df$d2< -.5))
  xyact=union(xact,yact)
  att_df$act=rep("inner",nrow(att_df))
  att_df$act[xyact]="outer"
  
  if(what[1]==TRUE && what[2]==FALSE ){
    a=ggplot(data=obs_df,aes(x=d1,y=d2))+coord_cartesian(xlim = c(xallmin,xallmax), ylim = c(yallmin,yallmax))
    a=a+geom_point(aes(x=d1,y=d2,colour=gr,shape=gr,alpha=.4),size=1)+theme_bw()
    #do not show obs labels if more than 30
    if (dim(x$obscoord)[1] < 30) {
      a=a+geom_text_repel(data=obs_df,aes(label=olab))
    }
    a=a+theme(legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank())+xlab("")+ylab("")
    if (length(x$size) != 1)
    {
      a=a+geom_point(data=group_df,aes(x=d1,y=d2,shape=gr))+theme(legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank())
      a=a+geom_text_repel(data=group_df,aes(label=glab))
    }
    a=a+geom_vline(xintercept=0)+geom_hline(yintercept=0)
    a=a+xlab(paste("Dim.",dims[1])) + ylab(paste("Dim.",dims[2]))  
    
    # if(disp==F){ggsave(filename = paste("K",deparse(K),"Map_units.pdf",sep=""),a,height=8 , width=8)
    #  }else{
    # out = a
    out$map=a
    #  print(a)
    #  }
    
  }
  if(what[1]==FALSE && what[2]==TRUE ){
    
    # xallmax=xattmax
    # xallmin=xattmin
    # yallmax=yattmax
    # yallmin=yattmin
    
    xallmax=1
    xallmin=-1
    yallmax=1
    yallmin=-1
    
    
    if(nrow(att_df)>=25){
      decr=(nrow(att_df)-25)*(1/250)
      mysize=5 * (1-decr)
      mysize=max(2,mysize)
    }else{mysize=5}
    
    a=ggplot(data=att_df,aes(x=d1,y=d2))+coord_cartesian(xlim = c(xallmin,xallmax), ylim = c(yallmin,yallmax))
    a=a+geom_point(alpha=.5,size=.25)+theme_bw()+xlab("")+ylab("")
    a=a+geom_text_repel(data=subset(att_df,act=="outer"),aes( label = attnam),size=mysize,segment.size = 0.1)
    a=a+geom_text(data=subset(att_df,act!="outer"),aes( label = attnam),size=mysize) #size=mysize*.8)
    a=a+geom_segment(data=att_df, aes(x=0,y=0,xend=d1,yend=d2),arrow=arrow(angle=15,unit(0.15, "inches")))
    a=a+annotate("path",x=0+1*cos(seq(0,2*pi,length.out=100)),
                 y=0+1*sin(seq(0,2*pi,length.out=100)))
    # a=a+geom_point(data=group_df,aes(x=d1,y=d2,shape=gr))+theme(legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank())
    # a=a+geom_text_repel(data=group_df,aes(label=glab))
    a=a+geom_vline(xintercept=0)+geom_hline(yintercept=0)
    a=a+xlab(paste("Dim.",dims[1])) + ylab(paste("Dim.",dims[2]))  
    #  if(disp==F){
    #    ggsave(filename = paste("K",deparse(K),"Map_attributes.pdf",sep=""),a,height=8 , width=8)
    #  }else{
    # out = a
    out$map=a
    #  print(a)
    #  }
  }
  if(what[1]==TRUE && what[2]==TRUE ){
    if(nrow(att_df)>=25){
      decr=(nrow(att_df)-25)*(1/250)
      mysize=5 * (1-decr)
      mysize=max(2,mysize)
    }else{mysize=5}
    
    a=ggplot(data=obs_df,aes(x=d1,y=d2))+coord_cartesian(xlim = c(xallmin,xallmax), ylim = c(yallmin,yallmax))
    a=a+geom_point(aes(x=d1,y=d2,shape=gr,alpha=.4),size=1)+theme_bw()#,colour=gr
    #do not show obs labels if more than 30
    if (dim(x$obscoord)[1] < 30) {
      a=a+geom_text_repel(data=obs_df,aes(label=olab))
    }
    
    a=a+theme(axis.text.x = element_blank(),axis.text.y = element_blank())+xlab("")+ylab("")
    if (length(x$size) != 1)
    {
      a=a+geom_point(data=group_df,aes(x=d1,y=d2,shape=gr))
      a=a+geom_text_repel(data=group_df,aes(label=glab))
    }
    a=a+geom_vline(xintercept=0)+geom_hline(yintercept=0)
    att_df$slp=att_df$d2/att_df$d1
    
    # arrow_df=data.frame(slp=att_df$slp)
    quad_check=sign(att_df[,1:2])
    marg_df=quad_check
    marg_mat=matrix(c(xallmin,yallmin,xallmax,yallmax),nrow=2)
    
    for(j in 1:2){
      neg_val=which(quad_check[j]<0)
      marg_df[neg_val,j]=marg_mat[j,1]
      marg_df[-neg_val,j]=marg_mat[j,2]
    }
    who_marg=apply(marg_df,1,function(x)which.min(abs(x)))
    
    arrow_df=data.frame(marg_df)
    for(i in 1:length(who_marg)){
      arrow_df$rd2[i]=arrow_df$d1[i]*(att_df$slp[i])
      arrow_df$rd1[i]=arrow_df$d2[i]*(1/att_df$slp[i])
    }
    sel_arrow_x=apply(arrow_df[,c(2,4)],1,function(x) which.min(abs(x)))
    
    myarrow_df=arrow_df[,1:2]
    for(i in 1:length(sel_arrow_x)){
      if(sel_arrow_x[i]==1){
        myarrow_df$d1[i]=arrow_df$d1[i]
        myarrow_df$d2[i]=arrow_df$rd2[i]
      }else{
        myarrow_df$d1[i]=arrow_df$rd1[i]
        myarrow_df$d2[i]=arrow_df$d2[i]
      }
    }
    
    myarrow_df$attnam=att_df$attnam
    
    a=a+geom_abline(data=att_df,aes(intercept=0,slope=slp,colour=attnam),alpha=.5)
    a=a+geom_segment(data=myarrow_df,aes(x=0,y=0,xend=d1,yend=d2,colour=attnam),alpha=.5,
                     arrow=arrow(length=unit(.15,"inches")))
    
    a=a+theme(legend.title=element_blank(),legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank())
    a=a+guides(shape=FALSE, alpha=FALSE)
    
    #do not show var labels if more than 50
    if (dim(x$attcoord)[1] < 50) {
      a=a+geom_text_repel(data=myarrow_df,aes(x=d1,y=d2,label=attnam))
    }
    
    a=a+xlab(paste("Dim.",dims[1])) + ylab(paste("Dim.",dims[2]))  
    out$map = a
    
  }
  print(a)
  if(cludesc==TRUE){
    #if the user gives custom attlabs, they won't be used
    #in cludesc=TRUE
    
    data <- as.data.frame(x$odata)
    data <- droplevels(data)
    numAct <- which(sapply(data, is.numeric))
    facAct <- which(!sapply(data, is.numeric))
    
    anynum <- any(numAct)
    anyfact <- any(facAct)
    

    if (anynum) {
      #FOR QUANTITATIVE
      QuantiAct <- as.matrix(data[, numAct, drop = FALSE])
      numobs = nrow(data)
      #standardize continuous
      QuantiAct <- t(t(QuantiAct) - as.vector(crossprod(rep(1,numobs)/numobs, as.matrix(QuantiAct))))
      QuantiAct <- t(t(QuantiAct)/sqrt(as.vector(crossprod(rep(1,numobs)/numobs, 
                                                           as.matrix(QuantiAct^2)))))
      X = QuantiAct
      cdsc = clu_means(X, x$cluster, center=FALSE, scale=FALSE)
      
      # out$map = a
      print(cdsc)
      out$parcoord = cdsc
    }
    
    if (anyfact){
      ### FOR CATEGORICAL
      
      csize = round((table(x$cluster)/sum(table(x$cluster)))*100,digits=1)
      cnames=paste("C",1:K,sep="")
      cnm=paste(cnames,": ",csize,"%",sep="")
      
      #QualiAct <-  tab.disjonctif(data[, facAct, drop = FALSE])
      
      lab1a=names(data[, facAct, drop = FALSE])
      lab1b=lapply(data[, facAct, drop = FALSE],function(z) levels(z))
      lab1=rep(lab1a,times=unlist(lapply(lab1b,length)))
      lab2=unlist(lab1b)
      qualilabs=paste(lab1,lab2,sep=".")
      
      attlabs=qualilabs #colnames(QualiAct)
      
      if (topstdres > length(attlabs)) {
        topstdres = length(attlabs)
      }
      ffew = topstdres 
      
      myminx = -10
      mymaxx = 10
      TopplotGroups=outOfIndependence(data[, facAct, drop = FALSE],x$cluster,attlabs,firstfew=ffew,textSize=4,segSize=4,minx=myminx,maxx=mymaxx)
      plotGroups=outOfIndependence(data[, facAct, drop = FALSE],x$cluster,nolabs=T,attlabs,fixmarg=F,textSize=1.5,segSize=1.5,minx=-2.5,maxx=2.5)
      
      for(jjj in 1:K){
        TopplotGroups$G[[jjj]]=TopplotGroups$G[[jjj]]+theme_bw()+ggtitle(cnm[jjj])
        
        if (subplot == TRUE) {
          out$stdres = TopplotGroups$G
          print(TopplotGroups$G[[jjj]])
          print(plotGroups$G[[jjj]], vp=viewport(.15, .18, .3, .35))
        }else{print(TopplotGroups$G[[jjj]])}
        # print(TopplotGroups$G[[jjj]])
      }
    }
  }
  invisible(out)
}

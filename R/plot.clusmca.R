plot.clusmca<-function(x, dims = c(1,2), what = c(TRUE,TRUE), cludesc = FALSE, topstdres = 20, objlabs = FALSE, attlabs = NULL, subplot = FALSE, ...){
  
  act = NULL
  attnam = NULL
  d1 = NULL
  d2 = NULL
  gr = NULL
  olab = NULL
  out=list()
  if (dim(data.frame(x$attcoord,stringsAsFactors = TRUE))[2] == 1) {
    stop('There is only one dimension. A 2D scatterplot cannot be produced.')
  } 
  
  dim1=dims[1]
  dim2=dims[2]
  K = max(x$cluster)
  
  dfAtt=data.frame(x1=x$attcoord[,1],x2=x$attcoord[,2])
  
  if (is.null(attlabs)) {
    lab1a=names(x$odata)
    lab1b=lapply(x$odata,function(z) levels(z))
    lab1=rep(lab1a,times=unlist(lapply(lab1b,length)))
    lab2=unlist(lab1b)
    attlabs=paste(lab1,lab2,sep=".")
  }
  
  
  #do not show obs labels if more than 30
  if (objlabs == TRUE) {
    obslabs = row.names(x$odata)
  } else
  {
    obslabs = paste("")
  }
  
  xallmax=max(max(x$attcoord[,dim1]),max(x$obscoord[,dim1]))
  xallmin=min(min(x$attcoord[,dim1]),min(x$obscoord[,dim1]))
  yallmax=max(max(x$attcoord[,dim2]),max(x$obscoord[,dim2]))
  yallmin=min(min(x$attcoord[,dim2]),min(x$obscoord[,dim2]))
  
  xallmax=max(max(x$attcoord[,dim1]),max(x$obscoord[,dim1]))
  xallmin=min(min(x$attcoord[,dim1]),min(x$obscoord[,dim1]))
  yallmax=max(max(x$attcoord[,dim2]),max(x$obscoord[,dim2]))
  yallmin=min(min(x$attcoord[,dim2]),min(x$obscoord[,dim2]))
  x_all_range=xallmax-xallmin
  y_all_range=yallmax-yallmin
  all_range=max(x_all_range,y_all_range)
  xallmax=xallmin+all_range
  yallmax=yallmin+all_range
  
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
  filt = 1*att_range
  att_df=data.frame(d1=x$attcoord[,dim1],d2=x$attcoord[,dim2],attnam=attlabs,stringsAsFactors = TRUE)
#  if(binary == TRUE){
#    pres=seq(from=2,to=nrow(att_df),by=2)
 #   print(pres)
#    att_df=att_df[pres,]
#  }
 
  xact=union(which(att_df$d1> filt),which(att_df$d1< -filt))
  yact=union(which(att_df$d2> filt), which(att_df$d2< -filt))
  xyact=union(xact,yact)
  att_df$act=rep("inner",nrow(att_df))
  att_df$act[xyact]="outer"
  
  glab=paste(rep("C",K),1:K,sep="")
  if (length(x$size) != 1)
  {
    group_df= data.frame(d1=x$centroid[,dim1],d2=x$centroid[,dim2],glab=glab,stringsAsFactors = TRUE)
  }
  obs_df=data.frame(d1=x$obscoord[,dim1],d2=x$obscoord[,dim2],gr=factor(x$cluster),olab=obslabs,stringsAsFactors = TRUE)
  
  if(what[1]==TRUE && what[2]==FALSE ){
    if (length(x$size) != 1)
    {
      names(group_df)[3] = "gr"
      levels(obs_df$gr) = levels(group_df$gr)
    }
    a=ggplot(data=obs_df,aes(x=d1,y=d2,colour=gr,shape=gr))+coord_cartesian(xlim=c(xallmin,xallmax),ylim=c(yallmin,yallmax))
    a=a+geom_point(aes(x=d1,y=d2,colour=gr,shape=gr,alpha=.4),size=1,na.rm = TRUE)+theme_bw()
    if (objlabs == TRUE) {
      a=a+geom_text_repel(data=obs_df,aes(label=olab))
    }
    
    a=a+theme(legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank())+xlab("")+ylab("")
    a=a+geom_vline(xintercept=0)+geom_hline(yintercept=0)
    
    if (length(x$size) != 1)
    {
      a=a+geom_point(data=group_df,colour="black",aes(x=d1,y=d2,shape=gr),na.rm=TRUE)+theme(legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank())
      
      a=a+geom_text_repel(data=group_df,colour="black",aes(label=gr))
    }
    a=a+xlab(paste("Dim.",dims[1])) + ylab(paste("Dim.",dims[2]))  
    #out = a
    out$map=a
    # print(a)
    
  }
  if(what[1]==FALSE && what[2]==TRUE ){
    
    xallmax=xattmax
    xallmin=xattmin
    yallmax=yattmax
    yallmin=yattmin
    
    if(nrow(att_df)>=25){
      decr=(nrow(att_df)-25)*(1/250)
      mysize=5 * (1-decr)
      mysize=max(2,mysize)
    }else{mysize=5}
    
    a=ggplot(data=att_df,aes(x=d1,y=d2))+coord_cartesian(xlim=c(xallmin,xallmax),ylim=c(yallmin,yallmax))
    a=a+geom_point(alpha=.5,size=.25,na.rm = TRUE)+theme_bw()+xlab("")+ylab("")
    a=a+geom_text_repel(data=subset(att_df,act=="outer"),aes( label = attnam),size=mysize,segment.size = 0.01)
    a=a+geom_text_repel(data=subset(att_df,act!="outer"),aes( label = attnam),size=mysize*.8,segment.size = 0.01)
    if (length(x$size) != 1)
    {  
      a=a+geom_point(data=group_df,aes(x=d1,y=d2,shape=glab),na.rm=TRUE)+theme(legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank())
      a=a+geom_text_repel(data=group_df,aes(label=glab))
    }
    a=a+geom_vline(xintercept=0)+geom_hline(yintercept=0)
    a=a+xlab(paste("Dim.",dims[1])) + ylab(paste("Dim.",dims[2]))  
    #out = a
    out$map=a
    # print(a)
  }
  if(what[1]==TRUE && what[2]==TRUE ){
    
    if (length(x$size) != 1)
    {
      names(group_df)[3] = "gr"
      levels(obs_df$gr) = levels(group_df$gr)
    }
    
    if(nrow(att_df)>=25){
      decr=(nrow(att_df)-25)*(1/250)
      mysize=5 * (1-decr)
      mysize=max(2,mysize)
    }else{mysize=5}
    
    a=ggplot(data=att_df,aes(x=d1,y=d2))+coord_cartesian(xlim=c(xallmin,xallmax),ylim=c(yallmin,yallmax))
    
    a=a+geom_point(data=obs_df,aes(x=d1,y=d2,colour=gr,shape=gr,alpha=.4),size=1,na.rm = TRUE)+theme_bw()
    if (objlabs == TRUE) {
      a=a+geom_text_repel(data=obs_df,aes(label=olab))
    }
    
    
    a=a+theme(legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank())+xlab("")+ylab("")
    a=a+geom_vline(xintercept=0)+geom_hline(yintercept=0)
    if (length(x$size) != 1)
    {
      a=a+geom_point(data=group_df,colour="black",aes(x=d1,y=d2,shape=gr),na.rm = TRUE)+theme(legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank())
      a=a+geom_text_repel(data=group_df,colour="black",aes(label=gr))
    }
    # 
    a = a + geom_point(data=att_df,aes(x=d1,y=d2),alpha=.5,size=.25,na.rm=TRUE) #+theme_bw()+xlab("")+ylab("")
    a=a+geom_text_repel(data=subset(att_df,act=="outer"),aes( label = attnam),size=mysize,segment.size = 0.1)
    a=a+geom_text_repel(data=subset(att_df,act!="outer"),aes( label = attnam),size=mysize*.8,segment.size = 0.01)
    a=a+geom_vline(xintercept=0)+geom_hline(yintercept=0)
    a=a+xlab(paste("Dim.",dims[1])) + ylab(paste("Dim.",dims[2]))  
    
  #######
    # a=a+geom_vline(xintercept=0)+geom_hline(yintercept=0)
    # att_df$slp=att_df$d2/att_df$d1
    # 
    # # arrow_df=data.frame(slp=att_df$slp)
    # quad_check=sign(att_df[,1:2])
    # marg_df=quad_check
    # marg_mat=matrix(c(xallmin,yallmin,xallmax,yallmax),nrow=2)
    # 
    # for(j in 1:2){
    #   neg_val=which(quad_check[j]<0)
    #   marg_df[neg_val,j]=marg_mat[j,1]
    #   marg_df[-neg_val,j]=marg_mat[j,2]
    # }
    # 
    # who_marg=apply(marg_df,1,function(x)which.min(abs(x)))
    # 
    # arrow_df=marg_df
    # for(i in 1:length(who_marg)){
    #   arrow_df$rd2[i]=arrow_df$d1[i]*(att_df$slp[i])
    #   arrow_df$rd1[i]=arrow_df$d2[i]*(1/att_df$slp[i])
    # }
    # 
    # sel_arrow_x=apply(arrow_df[,c(2,4)],1,function(x) which.min(abs(x)))
    # 
    # myarrow_df=arrow_df[,1:2]
    # for(i in 1:length(sel_arrow_x)){
    #   if(sel_arrow_x[i]==1){
    #     myarrow_df$d1[i]=arrow_df$d1[i]
    #     myarrow_df$d2[i]=arrow_df$rd2[i]
    #   }else{
    #     myarrow_df$d1[i]=arrow_df$rd1[i]
    #     myarrow_df$d2[i]=arrow_df$d2[i]
    #   }
    # }
    # 
    # myarrow_df$attnam=att_df$attnam
    # 
    # a=a+geom_abline(data=att_df,aes(intercept=0,slope=slp,colour=attnam),alpha=.5)
    # a=a+geom_segment(data=myarrow_df,aes(x=0,y=0,xend=d1,yend=d2,colour=attnam),alpha=.5,
    #                  arrow=arrow(length=unit(.15,"inches")))
    # 
    # a=a+theme(legend.title=element_blank(),legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank())
    # a=a+guides(shape=FALSE, alpha=FALSE)
    # 
    # #do not show var labels if more than 50
    # if (dim(x$attcoord)[1] < 50) {
    #   a=a+geom_text_repel(data=myarrow_df,aes(x=d1,y=d2,label=attnam))
    # }
    # 
    # a=a+xlab(paste("Dim.",dims[1])) + ylab(paste("Dim.",dims[2]))  
    # 
    
    
  ######
    out = list()
    out$map = a
    # print(a)
  }
  print(a)
  if(cludesc==TRUE){
    csize = round((table(x$cluster)/sum(table(x$cluster)))*100,digits=1)
    cnames=paste("C",1:K,sep="")
    cnm=paste(cnames,": ",csize,"%",sep="")
    
    if (topstdres > length(attlabs)) {
      topstdres = length(attlabs)
    }
    ffew = topstdres 
    
    if (subplot == TRUE) 
      TopplotGroups=outOfIndependence(x$odata,x$cluster,attlabs,firstfew=ffew,textSize=3.5,segSize=4,myleftmarg=5,myrightmarg=1)
    else
      TopplotGroups=outOfIndependence(x$odata,x$cluster,attlabs,firstfew=ffew,textSize=3.5,segSize=4,myleftmarg=0.5, myrightmarg=0.5)
    
    plotGroups=outOfIndependence(x$odata,x$cluster,nolabs=T,attlabs,fixmarg=F,textSize=1.5,segSize=1.5)#,myleftmarg=0.5, myrightmarg=0.5)
    
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
  
  invisible(out)
  
}
#}

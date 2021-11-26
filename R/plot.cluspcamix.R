plot.cluspcamix<-function(x, dims = c(1,2), cludesc = FALSE, topstdres = 20, objlabs = FALSE, attlabs = NULL, attcatlabs = NULL, subplot = FALSE, what = c(TRUE,TRUE), max.overlaps = 10, ...){
  d1 = NULL
  d2 = NULL
  gr = NULL
  olab = NULL
  act = NULL
  attnam = NULL
  slp = NULL
  out=list()
    
  data <- as.data.frame(x$odata,stringsAsFactors = TRUE)
  data <- droplevels(data)
  numAct <- which(sapply(data, is.numeric))
  facAct <- which(!sapply(data, is.numeric))
  
  if (dim(data.frame(x$attcoord,stringsAsFactors = TRUE))[2] == 1) {
    stop('There is only one dimension. A 2D scatterplot cannot be produced.')
  } 
  # out=list()
  dim1=dims[1]
  dim2=dims[2]
  K=max(x$cluster)
  
  if (is.null(attlabs)) {
    attlabs=rownames(x$attcoord)#colnames(x$odata) #row.names(x$attcoord)
  }
  
  if (is.null(attcatlabs)) {
    attcatlabs=rownames(x$attcatcoord)#colnames(x$odata) #row.names(x$attcoord)
  }
  
  if (objlabs == TRUE) {
    obslabs = row.names(x$odata)
  } else
  {
    obslabs = paste("")
  }
  
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
  
  xcatmax=max(max(x$attcatcoord[,dim1]),max(x$obscoord[,dim1]))
  xcatmin=min(min(x$attcatcoord[,dim1]),min(x$obscoord[,dim1]))
  ycatmax=max(max(x$attcatcoord[,dim2]),max(x$obscoord[,dim2]))
  ycatmin=min(min(x$attcatcoord[,dim2]),min(x$obscoord[,dim2]))
  x_cat_range=xcatmax-xcatmin
  y_cat_range=ycatmax-ycatmin
  cat_range=max(x_cat_range,y_cat_range)
  xcatmax=xcatmin+cat_range
  ycatmax=ycatmin+cat_range
  
  ######################################################
  att_df=data.frame(d1=x$attcoord[,dim1],d2=x$attcoord[,dim2],attnam=attlabs,stringsAsFactors = TRUE)
  att_catdf=data.frame(d1=x$attcatcoord[,dim1],d2=x$attcatcoord[,dim2],attnam=attcatlabs,stringsAsFactors = TRUE)
  
  glab=paste(rep("C",K),1:K,sep="")  
  if (length(x$size) != 1)
  {
    group_df= data.frame(d1=x$centroid[,dim1],d2=x$centroid[,dim2],glab=glab,gr=levels(factor(x$cluster)),stringsAsFactors = TRUE)
  }
  obs_df=data.frame(d1=x$obscoord[,dim1],d2=x$obscoord[,dim2],gr=factor(x$cluster),olab=obslabs,stringsAsFactors = TRUE)
  xact=union(which(att_df$d1> .5),which(att_df$d1< -.5))
  yact=union(which(att_df$d2>.5), which(att_df$d2< -.5))
  xyact=union(xact,yact)
  att_df$act=rep("inner",nrow(att_df))
  att_df$act[xyact]="outer"
  
  if(what[1]==TRUE && what[2]==FALSE ){
    a=ggplot(data=obs_df,aes(x=d1,y=d2))#+coord_cartesian(xlim = c(xallmin,xallmax), ylim = c(yallmin,yallmax))
    a=a+geom_point(aes(x=d1,y=d2,colour=gr,shape=gr,alpha=.4),size=1)+theme_bw()
    if (objlabs == TRUE) {
      a=a+geom_text_repel(data=obs_df,aes(label=olab), max.overlaps = max.overlaps)
    }
    a=a+theme(legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank())+xlab("")+ylab("")
    if (length(x$size) != 1)
    {
      a=a+geom_point(data=group_df,aes(x=d1,y=d2,shape=gr))+theme(legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank())
      a=a+geom_text_repel(data=group_df,aes(label=glab), max.overlaps = max.overlaps)
    }
    a=a+geom_vline(xintercept=0,color="grey")+geom_hline(yintercept=0,color="grey")
    a=a+xlab(paste("Dim.",dims[1])) + ylab(paste("Dim.",dims[2]))  
    
    # if(disp==F){ggsave(filename = paste("K",deparse(K),"Map_units.pdf",sep=""),a,height=8 , width=8)
    #  }else{
    # out = a
    out$map=a
    print(a)
    #  }
    
  }
  if(what[1]==FALSE && what[2]==TRUE ){
    
    a=ggplot(data=obs_df,aes(x=d1,y=d2))+theme_bw()#+coord_cartesian(xlim = c(xallmin,xallmax), ylim = c(yallmin,yallmax))
    # a=a+geom_point(aes(x=d1,y=d2,shape=gr,alpha=0),size=0)+theme_bw()#,colour=gr
    #do not show obs labels if more than 30
    # if (objlabs == TRUE) {
    #    a=a+geom_text_repel(data=obs_df,aes(label=olab), max.overlaps = max.overlaps)
    #  }
    
    a=a+theme(axis.text.x = element_blank(),axis.text.y = element_blank())+xlab("")+ylab("")
    if (length(x$size) != 1)
    {
      a=a+geom_point(data=group_df,aes(x=d1,y=d2,shape=gr))
      a=a+geom_text_repel(data=group_df,aes(label=glab), max.overlaps = max.overlaps)
    }
    a=a+geom_vline(xintercept=0,color="grey")+geom_hline(yintercept=0,color="grey")
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
    
    arrow_df=data.frame(marg_df,stringsAsFactors = TRUE)
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
    a=a+guides(shape="none", alpha="none")
    
    #do not show var labels if more than 50
    if (dim(x$attcoord)[1] < 50) {
      a=a+geom_text_repel(data=myarrow_df,aes(x=d1,y=d2,label=attnam),max.overlaps = max.overlaps)
    }
    
    a=a+xlab(paste("Dim.",dims[1])) + ylab(paste("Dim.",dims[2]))  
    out$map = a
    
    attquali_plot = att_catdf %>% ggplot(aes(x = d1, y = d2))+#geom_point(alpha=0.5)+
      #  theme(legend.title=element_blank(),legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank())+
      geom_text_repel(
        label=attcatlabs,max.overlaps = max.overlaps)+
      geom_vline(aes(xintercept=0),color="grey")+
      geom_hline(aes(yintercept=0),color="grey")+#coord_cartesian(xlim = c(xallmin,xallmax), ylim = c(yallmin,yallmax))+
      xlab(paste("Dim.",dims[1])) + ylab(paste("Dim.",dims[2]))  +
      theme_bw() +
      theme(legend.title=element_blank(),legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank())
    
    if (length(x$size) != 1)
    {
      attquali_plot=attquali_plot+geom_point(data=group_df,aes(x=d1,y=d2,shape=gr))
      attquali_plot=attquali_plot+geom_text_repel(data=group_df,aes(label=glab), max.overlaps = max.overlaps)
    }
    out$catmap <- attquali_plot
    out$map <- a
    print(attquali_plot)
    print(a)
    
    #  }
  }
  if(what[1]==TRUE && what[2]==TRUE ){
    a=ggplot(data=obs_df,aes(x=d1,y=d2))#+coord_cartesian(xlim = c(xCATmin,xallmax), ylim = c(yallmin,yallmax))
    a=a+geom_point(mapping=aes(x=d1,y=d2,shape=gr,alpha=.4),inherit.aes = FALSE,alpha=0.5,size=1)+theme_bw()
    if (objlabs == TRUE) {
      a=a+geom_text_repel(data=obs_df,aes(label=olab), max.overlaps = max.overlaps)
    }
    
    a=a+theme(axis.text.x = element_blank(),axis.text.y = element_blank())+xlab("")+ylab("")
    if (length(x$size) != 1)
    {
      a=a+geom_point(data=group_df,aes(x=d1,y=d2,shape=gr))
      a=a+geom_text_repel(data=group_df,aes(label=glab), max.overlaps = max.overlaps)
    }
    a=a+geom_vline(xintercept=0,color="grey")+geom_hline(yintercept=0,color="grey")
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
    
    arrow_df=data.frame(marg_df,stringsAsFactors = TRUE)
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
    a=a+guides(shape="none", alpha="none")
    
    #do not show var labels if more than 50
    if (dim(x$attcoord)[1] < 50) {
      a=a+geom_text_repel(data=myarrow_df,aes(x=d1,y=d2,label=attnam),max.overlaps = max.overlaps)
    }
    
    a=a+xlab(paste("Dim.",dims[1])) + ylab(paste("Dim.",dims[2]))  
    out$map = a
    #  scaling_factor = max(obs_df[,1]/att_df[,1])
    #  print(scaling_factor)
    #  print(head(att_catdf))
    #  att_catdf = scaling_factor*att_catdf[,1:2]
    
    
    attquali_plot = att_catdf %>% ggplot(aes(x = d1, y = d2))+
      # geom_text_repel(data=obs_df, mapping = aes(x = d1, y = d2, label=""), max.overlaps = 15, show.legend = FALSE,inherit.aes = F,alpha=1, size=4)+
      geom_text_repel(
        label=attcatlabs,max.overlaps = max.overlaps)+
      geom_vline(aes(xintercept=0),color="grey")+
      geom_hline(aes(yintercept=0),color="grey")+
      geom_point(data=obs_df,mapping = aes(x = d1, y = d2,shape=gr),inherit.aes = FALSE,alpha=0.5,size=1)+
      geom_text_repel(data=obs_df,aes(label=olab), max.overlaps = max.overlaps)+
      xlab(paste("Dim.",dims[1])) + ylab(paste("Dim.",dims[2]))  +
      theme_bw() +
      theme(legend.title=element_blank(),legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank())
    
    #   if (objlabs == TRUE) {
    #      a=a+geom_text_repel(data=obs_df,aes(label=olab), max.overlaps = max.overlaps)
    #    }
    
    if (length(x$size) != 1)
    {
      attquali_plot=attquali_plot+geom_point(data=group_df,aes(x=d1,y=d2,shape=gr))
      attquali_plot=attquali_plot+geom_text_repel(data=group_df,aes(label=glab), max.overlaps = max.overlaps)
    }
    
    out$catmap <- attquali_plot
    out$map <- a
    print(attquali_plot)
    print(a)
  }
  
  
  if(cludesc==TRUE){
    #if the user gives custom attlabs, they won't be used
    #in cludesc=TRUE
    
    
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
      
      out$map = a
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
      #  qualilabs = lab2
      attlabs=qualilabs #colnames(QualiAct)
      if (topstdres > length(attlabs)) {
        topstdres = length(attlabs)
      }
      ffew = topstdres 
      
      if (subplot == TRUE) 
        TopplotGroups=outOfIndependence(data[, facAct, drop = FALSE],Gvec=x$cluster,attlabs,firstfew=ffew,textSize=3.5,segSize=4,myleftmarg=5,myrightmarg=1)
      else
        TopplotGroups=outOfIndependence(data[, facAct, drop = FALSE],Gvec=x$cluster,attlabs,firstfew=ffew,textSize=3.5,segSize=4,myleftmarg=0.5, myrightmarg=0.5)
      
      
      # myminx = -10
      #  mymaxx = 10
      plotGroups=outOfIndependence(data[, facAct, drop = FALSE],Gvec=x$cluster,nolabs=TRUE,attlabs,fixmarg=FALSE,textSize=1.5,segSize=1.5)#,myleftmarg=0.5, myrightmarg=0.5)
      
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

# calculate gene signatures
calc_signatures <- function(x, gmtFile, method="mean", scale=F){
  
  geneset.obj<- GSA.read.gmt(gmtFile)
  genenames<-row.names(x)
  np=length(geneset.obj$genesets)
  
  if(scale){xs=t(scale(t(x),center=T,scale=T))}else{xs<-x}
  
  if(method!="gsa"){
    val=matrix(NA,nrow=np,ncol=ncol(x))
    for(i in 1:np){
      gene.set=which(genenames %in% geneset.obj$genesets[[i]])
      gene.set=gene.set[!is.na(gene.set)]
      if(length(gene.set)>1){
        if(method=="mean"){
          val[i,]=colSums(xs[gene.set,,drop=F])/length(gene.set)
        }
        if(method=="median"){
          val[i,]=t(apply(xs[gene.set,,drop=F],2,median))
        }
        if(method=="pca"){
          y<-prcomp(as.matrix(xs[gene.set,,drop=F]))
          val[i,]=y$rotation[,1]
        }
      } else if (length(gene.set) == 1) {
        val[i,] <- unlist(xs[gene.set,])
      }
    }
  }else{
    val<-GSA.make.features(gsaObj,xs,geneset.obj$genesets,genenames)
  }
  dimnames(val)<-list(geneset.obj$geneset.names,dimnames(x)[[2]])
  return(val)
}

overlapSets<-function(x,y){
  x<-x[dimnames(x)[[1]] %in% dimnames(y)[[1]],]
  y<-y[dimnames(y)[[1]] %in% dimnames(x)[[1]],]
  
  x<-x[sort.list(row.names(x)),]
  y<-y[sort.list(row.names(y)),]
  
  return(list(x=x,y=y))
}
assignDiffScore.dwd<-function(x,y){
  both<-overlapSets(x,y) 	# get the overlap of genes
  both$x<- apply(both$x,2,function(x){sign(x)*sqrt(x^2/sum(x^2))}) 	# set the distance to the origin to 1
  both$y<- apply(both$y,2,function(x){sign(x)*sqrt(x^2/sum(x^2))})	# set the distance to the origin to 1
  msproj<- apply(both$y,2,function(x,y){x%*%y},both$x[,1])	# project the samples on the MS-pL axis
  mlproj<- apply(both$y,2,function(x,y){x%*%y},both$x[,2])	# project the samples on the pL-mL axis
  diffScore<- mlproj - msproj 
  return( diffScore )	# return the point on the differentiation axis
}

GHI_RS <- function(edata) {
  getGene <- function(edata,gene){
    return(unlist(edata[rownames(edata)==gene,]))
  }
  gene <- c(2597,2990,60,7037,6175,2886,2064,596,5241,57758,2099,6790,4605,891,332,4288,4320,1515,968,2944,573)
  for(i in 1:length(gene)){
    assign(paste('X',gene[i],sep = ''),getGene(edata,gene[i]))
  }
  # normalized according to reference gene
  reference.Avg <- (X2597+X2990+X60+X7037+X6175)/5
  refNorm <- function(x,ref) {
    x <- x-ref
    x <- x-min(x)
    x <- x*15/max(x)
    return(x)
  }
  for(i in 1:length(gene)){
    assign(paste('X',gene[i],sep = ''),refNorm(get(paste('X',gene[i],sep = '')),reference.Avg))
  }
  
  GRB7_Group <- 0.9*X2886 + 0.1*X2064
  for(i in 1:ncol(edata))
  {
    if(GRB7_Group[i]<8)
    {GRB7_Group[i]<-8}
  }
  
  ER_Group<- (X596+1.2*X5241+X57758+0.8*X2099)/4
  
  Prolif_Group<- (X6790+X4605+X891+X332+X4288)/5
  
  for(i in 1:ncol(edata))
  {
    if(Prolif_Group[i]<6.5)
    {Prolif_Group[i]<-6.5}
  }
  
  Invasion_Group <- (X4320+X1515)/2
  
  CD68 <- X968
  GSTM1 <- X2944
  BAG1 <- X573
  
  RSU = 0.47*GRB7_Group - 0.34*ER_Group + 1.04*Prolif_Group + 0.10*Invasion_Group + 0.05*CD68 - 0.08*GSTM1 - 0.07*BAG1
  RS = 20*(RSU - 6.7)
  return(RS)
}

# test association between a signature and a gene using Fisher's exact test
fisherTest <- function(score,CN) {
  
  topQ <- quantile(score,0.75)
  table <- matrix(rep(0,4),nrow = 2,ncol = 2)
  module <- ifelse(score>=topQ,"high","low")
  temp <- data.frame(module=module,CN=CN)
  
  table[1,1] <- length(which(temp$module=="high" & temp$CN=="mut"))
  table[1,2] <- length(which(temp$module=="high" & temp$CN=="wt"))
  table[2,1] <- length(which(temp$module=="low" & temp$CN=="mut"))
  table[2,2] <- length(which(temp$module=="low" & temp$CN=="wt"))
  
  fish <- fisher.test(table,alternative = 'greater')
  return(c(fish$p.value,fish$estimate))
}

# test associations between a signature and all genes, both fisher and spearman correlation/lm
sigCNTest <- function(score,CN_score,CN_gain,CN_loss) {
  spearman_pos <- c()
  spearman_neg <- c()
  spearman_cor <- c()
  
  lm_pos <- c()
  lm_neg <- c()
  beta_coeff <- c()
  r.squared <- c()
  
  gain <- c()
  loss <- c()
  OR <- c()
  
  for(j in 1:nrow(CN_score)){
    CN <- unname(unlist(CN_score[j,]))
    
    # spearman rank correlation
    pos <- cor.test(score,CN,method = "spearman",alternative = "greater")
    neg <- cor.test(score,CN,method = "spearman",alternative = "less")
    spearman_pos <- c(spearman_pos,pos$p.value)
    spearman_neg <- c(spearman_neg,neg$p.value)
    spearman_cor <- c(spearman_cor,cor(score,CN,method = "spearman"))
    
    # linear model
    fit <- lm(score ~ CN + subtype_variable$basal + subtype_variable$her2 + subtype_variable$lumA + subtype_variable$lumB)
    sum <- summary(fit)
    beta <- sum$coefficients[,1][2]
    p <- sum$coefficients[,4][2]
    if(beta>0) {
      lm_pos <- c(lm_pos,p/2)
      lm_neg <- c(lm_neg,1-p/2)
    }
    else if(beta<0) {
      lm_pos <- c(lm_pos,1-p/2)
      lm_neg <- c(lm_neg,p/2)
    }
    beta_coeff <- c(beta_coeff,beta)
    r.squared <- c(r.squared,sum$r.squared)
  }
  
  # fisher's exact test
  gain <- apply(CN_gain,1,function(x){return(fisherTest(score,unlist(x)))})
  loss <- apply(CN_loss,1,function(x){return(fisherTest(score,unlist(x)))})
  
  p <- data.frame(spearman_pos = spearman_pos,spearman_neg = spearman_neg, CN_gain = gain[1,], CN_loss = loss[1,],lm_pos = lm_pos,lm_neg = lm_neg)
  return(p)
}

calc_segments<-function(x, gmtFile, method="mean", scale=F, gsaObj=NA){
  
  geneset.obj<- GSA.read.gmt(gmtFile)
  genenames<-row.names(x)
  np=length(geneset.obj$genesets)
  
  if(scale){xs=t(scale(t(x),center=T,scale=T))}else{xs<-x}
  
  if(method!="gsa"){
    val=matrix(0,nrow=np,ncol=ncol(x))
    for(i in 1:np){
      gene.set=which(genenames %in% geneset.obj$genesets[[i]])
      gene.set=gene.set[!is.na(gene.set)]
      if(length(gene.set)>1){
        if(method=="mean"){
          val[i,]=colSums(xs[gene.set,,drop=F])/length(gene.set)
        }
        if(method=="median"){
          val[i,]=t(apply(xs[gene.set,,drop=F],2,median))
        }
        if(method=="pca"){
          y<-prcomp(as.matrix(xs[gene.set,,drop=F]))
          val[i,]=y$rotation[,1]
        }
      } else if (length(gene.set) == 1) {
        val[i,] <- unlist(xs[gene.set,])
      }
    }
  }else{
    val<-GSA.make.features(gsaObj,xs,geneset.obj$genesets,genenames)
  }
  
  dimnames(val)<-list(geneset.obj$geneset.names,dimnames(x)[[2]])
  return(val)
}

medianCtr<-function(x){
  annAll <- dimnames(x)
  medians <- apply(x,1,median,na.rm=T)
  x <- t(scale(t(x),center=medians,scale=F))
  dimnames(x) <- annAll
  return(x)
}

standardize<-function(x){
  annAll<-dimnames(x)
  x<-scale(x)
  dimnames(x)<-annAll
  return(x)
}

exp_wrap <- function(edata) {
  keep <- c('29126','1493','5133') # for PD1, PDL1 and CTLA4
  edata70 <- edata[rowSums(edata<=2)<(0.3*ncol(edata)),]
  for(i in 1:3) {
    if(!(keep[i] %in% rownames(edata70))) {
      edata70 <- rbind(edata70,edata[keep[i],])
    }
  }
  edata70[edata70<=2] <- 0  # missing data marked with 0
  edata70log2 <- log(edata70,base = 2)
  edata70log2[edata70log2=="-Inf"] <- 0
  
  exp <- medianCtr(edata70log2)
  exp <- standardize(exp)
  return(exp)
}

caret_wrap <- function(trainX,trainY,testX,testY,bi=T) {
  if(!bi) {
    # set cross validation resampling method
    train_control <- trainControl(method = 'LGOCV',number = 200,classProbs = F)
    
    # set alpha and lambda grid
    alpha <- seq(0.1,0.9,by=0.1)
    lambda <- list()
    for(i in 1:9) {
      init <- glmnet(trainX,trainY,alpha = alpha[i])
      lambda[[i]] <- init$lambda
    }
    lambda_min <- min(unlist(lapply(lambda,min)))
    lambda_max <- max(unlist(lapply(lambda,max)))
    
    tune_grid = expand.grid(alpha = seq(0.1,0.9,by=0.1), lambda = seq(lambda_min,lambda_max,length.out = 100))
    
    # train model
    glmnet_obj <- train(trainX, trainY, method = "glmnet", metric = "RMSE",
                        trControl = train_control,tuneGrid = tune_grid)
    }  else {
    trainY2 <- ifelse(trainY==1,'pos','neg')
    
    # set cross validation resampling method
    train_control <- trainControl(method = 'LGOCV',number = 200,classProbs = T)
    
    alpha <- seq(0.1,0.9,by=0.1)
    lambda <- list()
    for(i in 1:9) {
      init <- glmnet(trainX,trainY,alpha = alpha[i],family = 'binomial')
      lambda[[i]] <- init$lambda
    }
    lambda_min <- min(unlist(lapply(lambda,min)))
    lambda_max <- max(unlist(lapply(lambda,max)))
    
    tune_grid = expand.grid(alpha = seq(0.1,0.9,by=0.1), lambda = seq(lambda_min,lambda_max,length.out = 100))
    
    # train model
    glmnet_obj <- train(trainX, trainY2, method = "glmnet", metric = "Accuracy",
                        trControl = train_control,tuneGrid = tune_grid)
    }
  return(glmnet_obj)
}

plot_ROC <- function(perf1,perf2,a1,a2,main) {
  tiff(paste0(main,'.tiff'),width = 1.5,height = 1.5,units = 'in',res = 300)
  par(mai = c(0.2,0.2,0.05,0.05),cex.axis = 0.3)
  plot(perf1,col = 'red',lwd = 1.2,xlab = "",ylab = "",box.lwd=0.8,
       xaxis.xaxt = 'n',yaxis.yaxt = 'n')
  plot(perf2,col = 'darkblue',lwd = 1.2,add = T)
  abline(a=0,b=1,lwd = 0.8)
  axis(1,tck=(-0.02),lwd = 0.8)
  axis(2,tck=(-0.02),lwd = 0.8)
  mtext(side = 1,text = seq(0,1,by=0.2),at = seq(0,1,by=0.2),cex = 0.5,line = (-0.2))
  mtext(side = 2,text = seq(0,1,by=0.2),at = seq(0,1,by=0.2),cex = 0.5,line = 0.1)
  
  legend('bottomright',legend = c(paste0('AUC = ',a1),paste0('AUC = ',a2)), lty = c(1,1),lwd = c(1,1) ,
         col = c('red','darkblue'),cex = 0.6,bty = 'n')
  dev.off()
}

x_y <- function(beta,segment_anno) {
  x <- c()
  y <- c()
  for(i in 1:nrow(segment_anno)) {
    this_seg <- rownames(segment_anno)[i]
    this_start <- segment_anno[i,2]
    this_end <- segment_anno[i,3]
    this_x <- seq(this_start,this_end,by = 1)
    this_y <- rep(beta[this_seg],length(this_x))
    x <- c(x,this_x)
    y <- c(y,this_y)
  }
  return(data.frame(x=x,y=y))
}

# plot model feature for single signature
plot_seg_ss <- function(beta,main) {
  total_gene <- 24776
  vertical <- vertical_lines[-c(1,24)]
  text_pos <- vertical[1]/2
  for(i in 2:length(vertical)){
    thispos <- vertical[i] - (vertical[i] - vertical[i-1])/2
    text_pos <- rbind(text_pos,thispos)
  }
  text_pos <- rbind(text_pos,total_gene-(total_gene-vertical[22])/2)
  
  min_y <- min(beta)
  max_y <- max(beta)
  # beta whole arm
  index <- grep('wholearm',names(beta))
  beta_wholearm <- beta[index]
  beta <- beta[-index]
  
  pos_beta <- beta[beta>0]
  neg_beta <- beta[beta<0]
  
  pos_seg <- names(pos_beta)
  neg_seg <- names(neg_beta)
  
  pos_anno <- segment_anno[pos_seg,]
  neg_anno <- segment_anno[neg_seg,]
  
  # pos regions excluding whole arms
  pos_coor <- x_y(beta,pos_anno)
  neg_coor <- x_y(beta,neg_anno)
  
  x <- c(pos_coor$x,neg_coor$x)
  y <- c(pos_coor$y,neg_coor$y)
  color <- c(rep('red',nrow(pos_coor)),rep('darkblue',nrow(neg_coor)))
  
  # whole arm beta
  beta_pos_wholearm <- beta_wholearm[beta_wholearm>0]
  beta_neg_wholearm <- beta_wholearm[beta_wholearm<0]
  
  if( length(beta_pos_wholearm) == 0) {
    pos_wholearm_coor <- data.frame(x=0,y=0)
  } else {
    pos_wholearm_anno <- segment_anno[names(beta_pos_wholearm),]
    pos_wholearm_coor <- x_y(beta_wholearm,pos_wholearm_anno)
  }
  
  if(length(beta_neg_wholearm)==0) {
    neg_wholearm_coor <- data.frame(x=0,y=0)
  } else {
    neg_wholearm_anno <- segment_anno[names(beta_neg_wholearm),]
    neg_wholearm_coor <- x_y(beta_wholearm,neg_wholearm_anno)
  }
  x_wholearm <- c(pos_wholearm_coor$x,neg_wholearm_coor$x)
  y_wholearm <- c(pos_wholearm_coor$y,neg_wholearm_coor$y)
  color_wholearm <- c(rep('pink',nrow(pos_wholearm_coor)),rep('lightblue',nrow(neg_wholearm_coor)))
  par(cex.axis = 2)
  png(filename = paste(main,'.png',sep = ''),width = 28,height = 5,res = 72,units = 'in')
  par(cex.axis = 2,cex.lab = 2.5,mai=c(0.6,1.5,0.6,0.5))
  plot(x_wholearm,y_wholearm,type = 'h',col = color_wholearm,xlim=c(1,24776),ylim=c(min_y,max_y),lwd = 1,xaxs="i",xaxt = 'n',ylab = '',xlab = "")
  abline(v=vertical_lines,lwd = 1)
  abline(h = 0,lwd = 1)
  lines(x,y,xaxt = "n",yaxt = "n",col = color,type = 'h')
  mtext(c(1:20,"",22,'x'),side = 1,at = text_pos,line = 1.5,cex = 2.5)
  mtext('weight',side = 2,at = 0,line = 4,cex = 2.5)
  dev.off()
}

P_Plot <- function(p,main,vertical_lines = v,y1 = NULL,y2 = NULL,label_up = NULL,label_down = NULL) {
  pos <- p[,1]
  neg <- p[,2]
  gain <- p[,3]
  loss <- p[,4]
  total_gene <- nrow(p)
  
  # remove non-significant lines
  pos_index <- which(!(pos>2 & gain>2))
  neg_index <- which(!(neg>2 & loss>2))
  pos[pos_index] <- 0
  gain[pos_index] <- 0
  neg[neg_index] <- 0
  loss[neg_index] <- 0
  
  m1 <- (max(pos,gain))
  m2 <- (max(neg,loss))
  
  # reformat log_p values for better plot
  index <- which(gain>pos)
  pos2 <- rep(0,total_gene)
  pos2[index] <- pos[index]
  
  index <- which(loss > neg)
  neg2 <- rep(0,total_gene)
  neg2[index] <- neg[index]
  
  if (is.null(y1) & is.null(y2)) {
    # init plot to get suitable axis labels
    print(paste('m1 =',m1,'; m2 =',m2))
  } else {
    vertical <- vertical_lines[-c(1,24)]
    text_pos <- vertical[1]/2
    for(i in 2:length(vertical)){
      thispos <- vertical[i] - (vertical[i] - vertical[i-1])/2
      text_pos <- rbind(text_pos,thispos)
    }
    text_pos <- rbind(text_pos,total_gene-(total_gene-vertical[22])/2)
    
    png(filename = paste(main,'.png',sep = ''),width = 28,height = 5,res = 72,units = 'in')
    par(cex.axis = 2,font.lab = 2)
    #layout(matrix(c(rep(1,2),2),ncol = 1))
    
    layout(matrix(c(1,2),ncol = 1),heights = c(y1+y2/6,y2+y1/6))
    
    par(mai=c(0,1.5,0.6,0.5),lwd = 2)
    plot(c(1:total_gene),pos,type = "h", xaxt = "n",yaxt = "n",lwd=2,xlim = c(1,24776),xaxs="i",yaxs='i',bty = '7',
         ylim = c(0,y1),col = "red",ylab = NA,xlab = NA)
    lines(c(1:total_gene),gain,xaxt = "n",yaxt = "n",col = "orange",lwd=2,type = 'h')
    lines(c(1:total_gene),pos2,xaxt = "n",yaxt = "n",col = "red",lwd=2,type = 'h')
    axis(2,labels = label_up,at = label_up,las = 1)
    abline(v = vertical_lines,lwd=1)
    abline(h = 2,lwd = 1,lty = 2)
    
    par(mai = c(0.6,1.5,0,0.5),lwd = 2)
    plot(c(1:total_gene),neg,type = "h", xaxt = "n",yaxt = "n",lwd=2,xlim = c(1,24776),xaxs="i",yaxs='i',bty = 'u',
         ylim = c(y2,0),col = "darkblue",ylab = NA,xlab = NA)
    lines(c(1:total_gene),loss,xaxt = "n",yaxt = "n",col = "cornflowerblue",lwd=2,type = 'h')
    lines(c(1:total_gene),neg2,xaxt = "n",yaxt = "n",col = "darkblue",lwd=2,type = 'h')
    abline(v = vertical_lines,lwd=1)
    abline(h = 2,lwd = 1,lty = 2)
    axis(2,labels = label_down,at = label_down,las = 1)
    
    mtext(c(1:20,"",22,'x'),side = 1,at = text_pos,line = 1.5,cex = 2.5)
    #mtext(0,side = 2,at = 0,line = 1,cex = 2)
    mtext(expression("-Log"[10]*'q'),side = 2,at = 0,line = 4,cex = 2.5)
    
    legend("bottom",fill = c("red","blue","orange","cornflowerblue"),
           legend = c(expression("-log"[10]*'q'[REGRESSION]*'(positive correlation)'),
                      expression("-log"[10]*'q'[REGRESSION]*'(negative correlation)'),
                      expression("-log"[10]*'q'[FISHER]*'(gain)'),
                      expression("-log"[10]*'q'[FISHER]*'(loss)')),bty = 'n',
           cex = 1, ncol = 2 )
    dev.off()
    
  }
}

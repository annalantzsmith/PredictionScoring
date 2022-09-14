#remove(list=ls())
#library(caret)
library(PRROC)
library(rlang)

## this function appears to have disappeared from the
##  loo 2.0.1 package;  was included in 2.0.0
kfold_split_balanced <- function(K = 10, x = NULL) {
  stopifnot(
    !is.null(x),
    length(unique(x)) == 2,
    K == as.integer(K),
    length(K) == 1,
    K > 1,
    K <= length(x)
  )
  x <- as.integer(as.factor(x)) - 1
  x0_id <- which(x == 0)
  x1_id <- which(x == 1)
  bins <- rep(NA, length(x))
  bins[x0_id] <- kfold_split_random(K = K, N = length(x0_id))
  bins[x1_id] <- kfold_split_random(K = K, N = length(x1_id))
  return(bins)
}

kfold_split_random <- function(K = 10, N = NULL) {
  stopifnot(
    !is.null(N),
    K == as.integer(K),
    N == as.integer(N),
    length(K) == 1,
    length(N) == 1,
    K > 1,
    K <= N
  )
  perm <- sample.int(N)
  idx <- ceiling(seq(from = 1, to = N, length.out = K + 1))
  bins <- .bincode(perm, breaks = idx, right = FALSE, include.lowest = TRUE)
  return(bins)
}

predscore <- function(x, ## model fitted to dataset 1
                      newdata = NULL, ## dataset 2 (for validation)
                      K = NULL, Ksize = NULL, #folds = NULL,
                      n.predict = NULL, saveModels=FALSE,
                      wtd = FALSE) {
  
  ## doesn't like if the formula in x contains unknown variables
  x$call$formula <- x$formula
  
  ## get dataset, number of rows
  olddata <- kfold_and_reloo_data(x)
  N <- nrow(olddata)
  
  ## get number of partitions (determined by olddata)
  if (is.null(K) & is.null(Ksize)){ K <- N
  } else if (is.null(Ksize)){ Ksize <- floor(N/K)
  } else { K <- floor(N/Ksize) }
  
  ## other set up
  if (is.null(n.predict)){ n.predict <- min(1000,N) }
  yname <- deparse(attr(terms(x),"variables")[[2]])
  isBin <- x$family$family=="binomial"
  
  ## sample indices for partitions
  if (isBin){
    folds.cv <- kfold_split_balanced(K = K, x = olddata[,yname])
    if (!is.null(newdata)){
      if (wtd) newdata <- getWtdSample( newdata, olddata, yname )
      folds.val <- kfold_split_balanced(K = K, x = newdata[,yname])
    } else { folds.val <- NULL }
  } else {
    folds.cv <- kfold_split_random(K = K,N = nrow(olddata))
    if (!is.null(newdata)){ folds.val <- kfold_split_random(K = K, N = nrow(newdata))
    } else { folds.val <- NULL }}
  
  ## summarise partitions
  sum1 <- data.frame( y = olddata[,yname],
                      k = as.factor(folds.cv) )
  sum2 <- data.frame( y = newdata[,yname],
                      k = as.factor(folds.val) )
  print(sum1 %>% group_by(k) %>% summarise(p_k = mean(y), n = n() ))
  print(sum2 %>% group_by(k) %>% summarise(p_k = mean(y), n = n() ))

  ## calculate scores
  predStuff <- lapply(1:K, calcScores, olddata=olddata, newdata=newdata,
                      x=x, folds.cv=folds.cv, folds.val=folds.val,
                      n.predict=n.predict, isBin=isBin, yname=yname, K=K,
                      saveModels=saveModels) #, wtd=wtd)
  results <- list( cv=cleanUp(predStuff,"cv",K=K,folds=folds.cv,isBin=isBin) )
  if (!is.null(newdata)){ results$val <- cleanUp(predStuff,"val",K=K,folds=folds.val,isBin=isBin)}
  if (saveModels){
    flats <- flatten(predStuff)
    results$fit_k <- flats[names(flats)=="fit_k"] }
  results$newdata <- newdata
  
  return(results) 
}

getWtdSample <- function( newdata, olddata, yname){
    p1 <- mean(olddata[,yname])
    p2 <- mean(newdata[,yname])
    
    newdata_y1_id <- which(newdata[,yname]==1)
    newdata_y0_id <- which(newdata[,yname]==0)
    
    if (p1 > p2){
      n_add <- (p1*nrow(newdata) - length(newdata_y1_id)) / ( 1 - p1 )
      newdata <- rbind(newdata, newdata[sample(newdata_y1_id,n_add,replace=TRUE),])
    } else {
      n_add <- (length(newdata_y1_id) - p1*nrow(newdata)) / (p1)
      newdata <- rbind(newdata, newdata[sample(newdata_y0_id,n_add,replace=TRUE),])
    }
    return(newdata)
}

calcScores <- function( k, olddata, newdata, x, folds.cv, folds.val, n.predict, isBin, yname, K, saveModels){ #, wtd ){
  omitted.cv <- which(folds.cv == k)
  if (!is.null(newdata)){ omitted.val <- which(folds.val == k) }
  
  #if (!is.null(newdata)){ if (wtd){
  #    Ksize <- floor(nrow(olddata/K))
  #    n.ones.val <- ceiling(mean(newdata[,yname])*Ksize)
  #    omitted.val <- c( sample(which(newdata[,yname]==1), n.ones.val),
  #                      sample(which(newdata[,yname]==0), Ksize - n.ones.val) )
  #  } else{ omitted.val <- which(folds.val == k) }}
  
  print(paste("k =",k,", hat{p}_k =", mean(olddata[omitted.cv,yname])))
  
  message("Fitting model ", k, " out of ", K)
  fit_k_call <- update(
    object = x,
    data = olddata[-omitted.cv,, drop=FALSE],
    subset = rep(TRUE, nrow(olddata) - length(omitted.cv)),
    weights = NULL,
    refresh = 0,
    open_progress = FALSE,
    evaluate = FALSE
  )
  if (!is.null(getCall(x)$offset)) {
    fit_k_call$offset <- x$offset[-omitted.cv]
  }
  fit_k_call$subset <- eval(fit_k_call$subset)
  fit_k_call$data <- eval(fit_k_call$data)
  ## fit the model
  capture.output(
    fit_k <- eval(fit_k_call)
  )
  
  ## cross-validation
  pred_k.cv <- getPred(fit_k,olddata[omitted.cv,, drop=FALSE],n.predict)
  predScores.cv <- qhat(olddata[omitted.cv,, drop=FALSE][,yname], pred_k.cv, isBin)
  
  ## validation
  if (!is.null(newdata)){ 
    pred_k.val <- getPred(fit_k,newdata[omitted.val,, drop=FALSE],n.predict)
    predScores.val <- qhat(newdata[omitted.val,, drop=FALSE][,yname], pred_k.val, isBin)
  } else { predScores.val <- NULL }
  
  results <- list(cv=predScores.cv, val=predScores.val)
  if (saveModels){ results$fit_k <- fit_k}
  return( results )
}

predscoreVal <- function(fit_k, ## models fitted to subsets of dataset 1
                      newdata, K, n.predict=NULL){ ## dataset 2 (for validation)
  
  ## other set up
  if (is.null(n.predict)){ n.predict <- min(1000,nrow(newdata)) }
  yname <- deparse(attr(terms(fit_k[[1]]),"variables")[[2]])
  isBin <- fit_k[[1]]$family$family=="binomial"
  
  ## sample indices for partitions
  if (K > 1){
    if (isBin){ folds.val <- kfold_split_balanced(K = K, x = newdata[,yname])
    } else { folds.val <- kfold_split_random(K = K,N = nrow(newdata)) }
  } else { folds.val <- rep(1,nrow(newdata)) }
  
  ## calculate scores
  predStuff <- lapply(1:K, calcScoresVal, newdata=newdata,
                      fit_ks=fit_k, folds.val=folds.val,
                      n.predict=n.predict, isBin=isBin, yname=yname, K=K)
  results <- list( val=cleanUp(predStuff,"val",K=K,folds=folds.val,isBin=isBin) )
  
  return(results) 
}

calcScoresVal <- function( k, newdata, fit_ks, folds.val, n.predict, isBin, yname, K){
  omitted.val <- which(folds.val == k)
  fit_k <- fit_ks[[1]]
  pred_k.val <- getPred(fit_k,newdata[omitted.val,, drop=FALSE],n.predict)
  predScores.val <- qhat(newdata[omitted.val,, drop=FALSE][,yname], pred_k.val, isBin)
  results <- list(val=predScores.val)
  return( results )
}

cleanUp <- function(predStuff,mode,K,folds,isBin){
  firstFlats <- unlist(predStuff,recursive=FALSE)
  if (mode=="cv"){ flats <- flatten(firstFlats[names(firstFlats)=="cv"])
  } else if (mode=="val"){ flats <- flatten(firstFlats[names(firstFlats)=="val"]) }
  results <- data.frame( set=1:K,
                         q=unlist(flats[names(flats)=="q"]))
  if (isBin){
    results$q.prec2 <- unlist(flats[names(flats)=="q.prec2"])
    results$q.roc <- unlist(flats[names(flats)=="q.roc"])
    precCurves <- flats[names(flats)=="prec"]
    rocCurves <- flats[names(flats)=="roc"]
  } else { precCurves <- NULL; rocCurves <- NULL }
  idxSets <- lapply(1:K, function(k) which(folds == k))
  return(list(results=results,rocCurves=rocCurves,precCurves=precCurves,idxSets=idxSets))
}

getPred <- function(model,dataset,n=NA){
  #formulaWords <- unlist(strsplit(as.character(model$formula)," "))
  #modelTerms <- unique(formulaWords[!formulaWords %in% (c("~","+","*","1"))])
  preddata <- get_all_vars(model$formula,dataset)
  posterior_predict(model,
                    newdata=preddata,
                    draws=n)
}

qhat <- function( obs, dist, isBin ){
  if (!isBin){
    if (length(obs)>1) qhat.vector(obs,dist,isBin)
    else qhat.scalar(obs,dist,isBin)
  } else {
    aucStat(obs,dist)
  }
}

qhat.scalar <- function( obs, dist, isBin ){
  if (isBin) return( list(q=sum(obs==dist)/length(dist)) )
  else return( list(q=sum(obs>dist)/length(dist)) )
}

qhat.vector <- function( obs, dist, isBin ){
  if (length(obs)!=ncol(dist)) return("error: length(obs) must equal ncol(dist)")
  else return( list(q=sapply(1:length(obs),
                             function(x){ qhat.scalar(obs[x],dist[,x],isBin)$q })) )
}

ksVec <- function( paramMat1, paramMat2, type="statistic" ){
  sapply( 1:ncol(paramMat1), function(i){
    as.numeric(ks.test(paramMat1[,i],paramMat2[,i])[type]) })
}

aucStat <- function( truey, phat ){
  if (!is.list(phat)){ phat <- colMeans(phat) }
  fg <- phat[truey == 1]
  bg <- phat[truey == 0]
  roc <- roc.curve(scores.class0=fg, scores.class1=bg, curve=TRUE)
  prec <- pr.curve(scores.class0=fg, scores.class1=bg, curve=TRUE)
  #threshes <- prediction(phat,truey)
  precCurve <- data.frame( rec=prec$curve[,1],
                           prec=prec$curve[,2],
                           auc=rep(prec$auc.integral,nrow(prec$curve)),
                           auc2=rep(prec$auc.davis.goadrich,nrow(prec$curve)) )
  #tpr=performance(threshes,measure="tpr")@y.values[[1]],
  #fpr=performance(threshes,measure="fpr")@y.values[[1]],
  #prec=performance(threshes,measure="prec")@y.values[[1]],
  #rec=performance(threshes,measure="rec")@y.values[[1]])
  rocCurve <- data.frame( fpr=roc$curve[,1],
                          tpr=roc$curve[,2],
                          auc=rep(roc$auc,nrow(roc$curve)) )
  #auc <- performance(threshes,measure="auc")
  #precrecf <- performance(threshes,measure=)
  return(list( q=prec$auc.integral, #q=unlist(auc@y.values),
               q.prec2=prec$auc.davis.goadrich,
               q.roc=roc$auc,
               prec=precCurve, roc=rocCurve))
  ## want to also return a list of the pr objects to make matching roc curves
}

kfold_and_reloo_data <- function(x) {
  d <- x[["data"]] # either data frame or environment
  
  sub <- getCall(x)[["subset"]]
  if (!is.null(sub)) {
    keep <- eval(substitute(sub), envir = d)
  }
  
  if (is.environment(d)) {
    d <- get_all_vars(formula(x), d) # now d is a data frame
  }
  if (!is.null(sub)) {
    d <- d[keep,, drop=FALSE]
  }
  
  #d <- na.omit(d)
  getNArows <- attr(na.omit(d[,names(attr(x$terms,"dataClasses"))]),
                   "na.action")[1]
  d <- d[!(1:nrow(d) %in% getNArows),]
  
  if (is_clogit(x)) {
    strata_var <- as.character(getCall(x)$strata)
    d[[strata_var]] <- model.weights(model.frame(x))
  }
  return(d)
}

is_clogit <- function(object) {
  is(object, "clogit")
}

plot.predscore <- function( q1, q2=NULL,
                            main1="Cross-Validation",
                            main2="Validation",
                            plotNorm=FALSE,
                            isBin=FALSE,
                            plotHists=TRUE,
                            plot=TRUE){
  if (isBin){ plot.predscore.isBin(q1,q2,plotHists,plot)
  } else { plot.predscore.notBin(q1,q2,main1,main2,plotNorm) } #,plotHists,plot) }
}

plot.predscore.isBin <- function(cv,val,
                                 plotHists=TRUE,
                                 plot=TRUE,
                                 colorlist=c( brewer.pal(7,"BrBG")[6],
                                              "orangered")){
  q1 <- cv$results$q
  q2 <- val$results$q
  rocs1 <- cv$rocCurves
  rocs2 <- val$rocCurves
  precs1 <- cv$precCurves
  precs2 <- val$precCurves
  
  ## Unlist Curves
  rocdf <- NULL
  for (i in 1:length(rocs1)){ rocs1[[i]]$type <- as.factor("cross-validation") }
  for (i in 1:length(rocs2)){ rocs2[[i]]$type <- as.factor("validation") }
  rocs <- c(rocs1,rocs2)
  for (i in 1:length(rocs)){
    rocs[[i]]$lines <- i
    rocdf <- rbind(rocdf,rocs[[i]])
  }
  precdf <- NULL
  for (i in 1:length(precs1)){ precs1[[i]]$type <- as.factor("cross-validation") }
  for (i in 1:length(precs2)){ precs2[[i]]$type <- as.factor("validation") }
  precs <- c(precs1,precs2)
  for (i in 1:length(precs)){
    precs[[i]]$lines <- i
    precdf <- rbind(precdf,precs[[i]])
  }
  
  ## Find bands on curves
  rocbands <- rocdf %>% group_by(type,fpr) %>% 
    summarise( ymax=max(tpr),
               ymin=min(tpr),
               ymean=mean(tpr),
               ymed=median(tpr))
  precbands <- precdf %>% group_by(type,rec) %>% 
    summarise( ymax=max(prec),
               ymin=min(prec),
               ymean=mean(prec),
               ymed=median(prec) )  
  
  ## Make Histogram dfs
  qsdf <- data.frame(q=c(q1,q2),
                     type=c( rep("cross-validation",length(q1)),
                             rep("validation",length(q2)) ))
  qsdf2 <- data.frame( q = unlist(lapply(1:max(rocdf$lines),function(i){
    rocdf$auc[rocdf$lines==i][1] })),
    type = unlist(lapply(1:max(rocdf$lines),function(i){
      rocdf$type[rocdf$lines==i][1] })) )
  
  plotList <- list()
  ## Histograms
  plotList[[3]] <- ggplot(data=qsdf) + theme_bw() + theme(legend.position="bottom",legend.text=element_text(size=9)) +
    geom_histogram(aes(x=q,fill=type,col=type),position="identity",bins=25,alpha=.25) +
    labs(x="AUC statistics (PR)",y="Frequency") + 
    scale_fill_manual(values=colorlist,aesthetics=c("fill","color"))
  plotList[[1]] <- ggplot(data=qsdf2) + theme_bw() + theme(legend.position="bottom",legend.text=element_text(size=9)) +
    geom_histogram(aes(x=q,fill=type,col=type),position="identity",alpha=.25) +
    labs(x="AUC statistics (ROC)",y="Frequency") + 
    scale_fill_manual(values=colorlist,aesthetics=c("fill","color"))
  ## Curves
  plotList[[4]] <- ggplot() + theme_bw() + theme(legend.position="bottom",legend.text=element_text(size=9)) +
    geom_line(data=precdf,aes(x=rec,y=prec,group=lines,col=type),alpha=.10) +
    geom_ribbon(data=precbands,aes(rec,ymin=ymin,ymax=ymax,fill=type),alpha=.25) +
    geom_line(data=precbands,aes(rec,ymed,col=type),lwd=1) +
    labs(x="Recall",y="Precision") + 
    scale_fill_manual(values=colorlist,aesthetics=c("fill","color"))
  plotList[[2]] <- ggplot() + theme_bw() + theme(legend.position="bottom",legend.text=element_text(size=9)) +
    geom_line(data=rocdf,aes(x=fpr,y=tpr,group=lines,col=type),alpha=.10) +
    geom_ribbon(data=rocbands,aes(fpr,ymin=ymin,ymax=ymax,fill=type),alpha=.25) +
    geom_line(data=rocbands,aes(fpr,ymed,col=type),lwd=1) +
    labs(x="False positive rate",y="True positive rate") + 
    scale_fill_manual(values=colorlist,aesthetics=c("fill","color"))

  ## Boxplots
  plotList[[6]] <- ggplot(data=qsdf) + theme_bw() + theme(legend.position="bottom",legend.text=element_text(size=9)) +
    geom_boxplot(aes(y=q,fill=type,col=type),alpha=.25) +
    labs(x="AUC statistics (PR)",y="Frequency") + 
    scale_fill_manual(values=colorlist,aesthetics=c("fill","color"))
  plotList[[5]] <- ggplot(data=qsdf2) + theme_bw() + theme(legend.position="bottom",legend.text=element_text(size=9)) +
    geom_boxplot(aes(y=q,fill=type,col=type),alpha=.25) +
    labs(x="AUC statistics (ROC)",y="Frequency") + 
    scale_fill_manual(values=colorlist,aesthetics=c("fill","color"))
  
  ## Boxplots
  plotList[[8]] <- ggplot(data=qsdf) + theme_bw() + theme(legend.position="bottom",legend.text=element_text(size=9)) +
    geom_violin(aes(x=type,y=q,fill=type,col=type),alpha=.25,draw_quantiles=c(.25,.48,.49,.51,.52,.75)) +
    geom_jitter(aes(x=type,y=q,col=type),width=.1) +
    labs(x="AUC statistics (PR)",y="Frequency") + ylim(range(c(qsdf$q,qsdf2$q))) +
    scale_fill_manual(values=colorlist,aesthetics=c("fill","color")) 
  plotList[[7]] <- ggplot(data=qsdf2) + theme_bw() + theme(legend.position="bottom",legend.text=element_text(size=9)) +
    geom_violin(aes(x=type,y=q,fill=type,col=type),alpha=.25,draw_quantiles=c(.25,.48,.49,.51,.52,.75)) +
    geom_jitter(aes(x=type,y=q,col=type),width=.1) +
    labs(x="AUC statistics (ROC)",y="Frequency") + ylim(range(c(qsdf$q,qsdf2$q))) +
    scale_fill_manual(values=colorlist,aesthetics=c("fill","color"))
  
  ## Histogram of differences
  diffs_pr <- q2-q1
  plotList[[10]] <- ggplot() + theme_classic() + theme(legend.position="bottom",legend.text=element_text(size=9)) +
    geom_histogram(aes(diffs_pr,stat(density)),binwidth=0.04,fill=colorlist[2],col=colorlist[2],alpha=.25,lwd=.1) + 
    geom_density(aes(diffs_pr),col=colorlist[2],lwd=0.75) + 
    geom_vline(xintercept=0,lty=2) + geom_vline(xintercept=mean(diffs_pr),col=colorlist[2],lwd=0.75) +
    labs(x="Difference in AUC statistics (PR)",y="Frequency")
  diffs_roc <- filter(qsdf2,type=="cross-validation")$q - filter(qsdf2,type=="validation")$q
  plotList[[9]] <- ggplot() + theme_classic() + theme(legend.position="bottom",legend.text=element_text(size=9)) +
    geom_histogram(aes(diffs_roc,stat(density)),binwidth=0.04,fill=colorlist[2],col=colorlist[2],alpha=.25,lwd=.1) + 
    geom_density(aes(diffs_roc),col=colorlist[2],lwd=0.75) + 
    geom_vline(xintercept=0,lty=2) + geom_vline(xintercept=mean(diffs_roc),col=colorlist[2],lwd=0.75) +
    labs(x="Difference in AUC statistics (ROC)",y="Frequency")
  
  if (!plotHists){ plotList <- plotList[c(2,4)] }
  if (plot){ do.call("grid.arrange",c(plotList,ncol=length(plotList)))
  } else { return(plotList) }
}

plot.predscore.notBin <- function(cv,val,main1,main2,plotNorm){
  
  q1 <- cv$results$q
  q2 <- val$results$q
  
  if (length(q2)>0){ par(mfrow=c(1,3))
  } else { par(mfrow=c(1,2)) }
  
  nComp <- length(q1)
  if (length(q2)>0){ nComp <- min(c(length(q1),length(q2))) }
  
  compDist <- runif(nComp)
  compName <- "Uniform"
  xlabs <- "quantiles"
  xlim1 <- xlim2 <- c(0,1)
  if (plotNorm){
    q1 <- qnorm(q1); q1[q1==-Inf] <- NA; q1[q1==Inf] <- NA
    if (length(q2)>0){ q2 <- qnorm(q2); q2[q2==-Inf] <- NA; q2[q2==Inf] <- NA } 
    compDist <- rnorm(nComp)
    compName <- "Normal"
    xlabs <- "transformed quantiles"
    xlim1 <- range(q1,na.rm=TRUE); xlim2 <- range(q2,na.rm=TRUE)
  }
  
  ## Raw Quantiles (Uniform)
  hist(q1,breaks=25,main=main1,xlim=xlim1,
       probability=TRUE,xlab=xlabs)
  if (length(q2)>0){ hist(q2,breaks=25,main=main2,xlim=xlim2,
                          probability=TRUE,xlab=xlabs) }
  
  plot(ecdf(q1),main="Empirical CDFs",xlab=xlabs,ylab="",
       do.points=FALSE,verticals=TRUE,col="green")
  if (length(q2)>0){ plot(ecdf(q2),do.points=FALSE,verticals=TRUE,
                          col="blue",add=TRUE) }
  plot(ecdf(compDist),do.points=FALSE,verticals=TRUE,
       col="red",add=TRUE)
  if (length(q2)>0){ legend("topleft",c(main1,main2,compName),
                            col=c("green","blue","red"), lty=1, cex=.5)
  } else { legend("topleft",c(main1,compName),
                  col=c("green","red"), lty=1, cex=.5) }
}

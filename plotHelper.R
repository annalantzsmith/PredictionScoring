## Helper functions for the Gallup data

## Function to prepare round-level data
prepRoundData <- function( gallupdata, buffer=.01 ){
  session_info = gallupdata %>%
    filter(round_num==1) %>%
    group_by(sessionnum) %>%
    summarise(
      num_player = n(),
      condition = unique(condition)[1]
    ) %>%
    arrange(condition)
  
  session_round_rate = gallupdata%>%
    group_by(sessionnum, 
             round_num)%>%
    summarise(
      rate_contr = mean(decision..0.D.1.C.)
    )
  session_round_rate = left_join(session_round_rate, 
                                 session_info,
                                 by="sessionnum")
  
  boxplot_info = session_round_rate %>%
    group_by(round_num,condition) %>%
    summarise(
      avg_rate = mean(rate_contr),
      Q1 = quantile(rate_contr,probs=.25),
      Q3 = quantile(rate_contr,probs=.75),
      IQR = Q3 - Q1,
      highbar = Q3 + 1.5*IQR,
      lowbar = Q1 - 1.5*IQR,
      gghighbar = max(min(highbar,max(rate_contr)) - buffer,Q3),
      gglowbar = min(max(lowbar,min(rate_contr)) + buffer,Q1)
    )
  
  return(list(session_info=session_info,
              boxplot_info=boxplot_info,
              session_round_rate=session_round_rate))
}

# plot contribution rates
plotContr <- function( gallupRounds, title="" ){
  session_round_rate <- gallupRounds$session_round_rate
  boxplot_info <- gallupRounds$boxplot_info
  g <- ggplot(session_round_rate, 
              aes(x=factor(round_num), 
                  y=rate_contr,
                  fill=condition)) +
    ggtitle(title) + xlab("Round") + ylab("Contribution rate") + labs(fill="Network\nCondition",col=NA) +
    geom_line(aes(col=condition,group=sessionnum),show.legend=FALSE) +
    geom_boxplot(alpha=.5,varwidth=TRUE,col="gray25",coef=0,outlier.color=NA) +
    geom_segment(data=boxplot_info,aes(x=factor(round_num),y=gglowbar,
                                       xend=factor(round_num),yend=Q1),col="gray25") +
    geom_segment(data=boxplot_info,aes(x=factor(round_num),y=Q3,
                                       xend=factor(round_num),yend=gghighbar),col="gray25") +
    facet_grid(.~condition) +
    theme(axis.text.x=element_text(angle=0, vjust=0.4, hjust=1)) + 
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black"),
          strip.text = element_text(size=12)) +
    ylim(c(0,1))
  return(g)
}

## calculate ROC and PR curves
getCurves <- function( model, pred_data, 
                       yname = "decision..0.D.1.C.",
                       old_data=NULL, wtd=FALSE){
  if (wtd) pred_data <- getWtdSample( pred_data, old_data, yname)
  preds <- predict( model, pred_data )
  roc <- roc.curve( preds[pred_data[,yname]==1],
                    preds[pred_data[,yname]==0], curve=TRUE)
  rocdf <- as.data.frame(roc$curve)
  rocdf$wtd <- wtd
  
  pr <- pr.curve( preds[pred_data[,yname]==1],
                  preds[pred_data[,yname]==0], curve=TRUE)
  prdf <- as.data.frame(pr$curve)
  prdf$wtd <- wtd
  return( list( preds=preds, roc=roc, pr=pr, data=pred_data,
                rocdf=rocdf, prdf=prdf ))
}

## plot model estimates (regression line) with data
plotModel <- function( sessiondata, coef, title="", buffer=.01 ){
  sessiondata <- sessiondata %>%
    mutate(isfluid = (condition=="Fluid"),
           phat = invlogit( coef[1] + coef[2]*isfluid + coef[3]*round_num +
                              coef[4]*isfluid*round_num) )
  
  ggplot(sessiondata, 
         aes(x=factor(round_num), 
             y=rate_contr,
             fill=isfluid)) +
    ggtitle(title) + xlab("Round") + ylab("Contribution rate") + labs(fill="Fluid\nCondition",col=NA) +
    geom_boxplot(alpha=.15,varwidth=TRUE,col="gray25") + 
    geom_line(aes(x=round_num,y=phat,col=isfluid), lwd=2) +
    theme(axis.text.x=element_text(angle=0, vjust=0.4, hjust=1)) + 
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black"),
          strip.text = element_text(size=12))
}